#!/usr/bin/env python3
"""NTRIP test client for nrtk-svr (spec ch.11/12 real-time checks).

Connects to a caster mountpoint (NTRIP v1 or v2), sends a GGA sentence
periodically, receives RTCM 3 and reports:
  - time from connection to the first MSM epoch (target <= 3 s)
  - age of corrections per epoch (receive time - MSM epoch time, p95 <= 1.5 s;
    in replay mode pass --replay-start so the replay clock can be modelled)
  - message type counts, station ids, CRC errors
Optionally writes the received bytes to a file (decode with RTKLIB convbin).

usage:
  ntrip_client.py --host 127.0.0.1 --port 2101 --mount VRS_MSM7 \\
                  --llh 24.78 121.15 250 [--seconds 60] [--out vrs.rtcm3] \\
                  [--clients 100] [--v2] [--user u --password p]
"""
import argparse, base64, datetime, json, math, socket, statistics, struct, sys, threading, time

GPS0 = datetime.datetime(1980, 1, 6)
LEAP = 18


def crc24q(data):
    crc = 0
    for b in data:
        crc ^= b << 16
        for _ in range(8):
            crc <<= 1
            if crc & 0x1000000:
                crc ^= 0x1864CFB
    return crc & 0xFFFFFF


def getbits(buf, pos, n):
    v = 0
    for i in range(pos, pos + n):
        v = (v << 1) | ((buf[i // 8] >> (7 - i % 8)) & 1)
    return v


def gga_sentence(lat, lon, h, t=None):
    t = t or datetime.datetime.utcnow()
    def dm(x, w):
        a = abs(x); d = int(a); m = (a - d) * 60
        return '%0*d%07.4f' % (w, d, m)
    body = 'GPGGA,%02d%02d%05.2f,%s,%s,%s,%s,1,12,1.0,%.3f,M,0.000,M,,' % (
        t.hour, t.minute, t.second + t.microsecond / 1e6, dm(lat, 2), 'N' if lat >= 0 else 'S',
        dm(lon, 3), 'E' if lon >= 0 else 'W', h)
    cs = 0
    for c in body:
        cs ^= ord(c)
    return '$%s*%02X\r\n' % (body, cs)


def http_status(host, port):
    """GET /status of nrtk-svr -> (json, wall time of the reply)"""
    s = socket.create_connection((host, port), timeout=5)
    s.sendall(('GET /status HTTP/1.1\r\nHost: %s\r\nConnection: close\r\n\r\n' % host).encode())
    buf = b''
    while True:
        d = s.recv(65536)
        if not d:
            break
        buf += d
    t = time.time()
    s.close()
    body = buf.split(b'\r\n\r\n', 1)[1]
    return json.loads(body.decode()), t


class Client(threading.Thread):
    def __init__(self, a, idx):
        super().__init__(daemon=True)
        self.a, self.idx = a, idx
        self.types, self.staids = {}, set()
        self.ages, self.first_msm, self.crc_err, self.nbytes = [], None, 0, 0
        self.error, self.closed_by_server = None, False
        self.raw = bytearray()
        self.timeline = []   # (s since connect, type, staid, tow)

    def gps_now(self):
        now = datetime.datetime.utcnow()
        sec = (now - GPS0).total_seconds() + LEAP
        if self.a.replay_start is not None:   # model the replay clock
            w0 = getattr(self.a, 'replay_wall0', None) or self.t_connect
            sec = self.a.replay_start + (time.time() - w0) * self.a.speed + self.a.replay_offset
        return sec % 604800.0

    def run(self):
        a = self.a
        try:
            s = socket.create_connection((a.host, a.port), timeout=10)
        except OSError as e:
            self.error = str(e)
            return
        self.t_connect = time.time()
        hdr = 'GET /%s HTTP/1.1\r\nHost: %s\r\nUser-Agent: NTRIP nrtk-test/1.0\r\n' % (a.mount, a.host)
        if a.v2:
            hdr += 'Ntrip-Version: Ntrip/2.0\r\n'
        if a.user:
            hdr += 'Authorization: Basic %s\r\n' % base64.b64encode(('%s:%s' % (a.user, a.password)).encode()).decode()
        s.sendall((hdr + '\r\n').encode())
        lat, lon, h = a.llh
        lat += a.spread * math.sin(self.idx * 2.39996) / 111.0
        lon += a.spread * math.cos(self.idx * 2.39996) / 111.0
        if not a.no_gga:
            s.sendall(gga_sentence(lat, lon, h).encode())
        last_gga = time.time()
        buf = b''
        header_done, chunked = False, False
        s.settimeout(0.5)
        frames = bytearray()
        while time.time() - self.t_connect < a.seconds:
            try:
                d = s.recv(65536)
            except socket.timeout:
                d = None
            except OSError as e:
                self.error = str(e)
                break
            if d == b'':
                self.closed_by_server = True
                break
            if d:
                buf += d
            if not header_done and buf.startswith(b'ICY 200 OK\r\n'):   # NTRIP v1: no blank line
                buf = buf[12:]
                if buf.startswith(b'\r\n'):
                    buf = buf[2:]
                header_done = True
            elif not header_done and b'\r\n\r\n' in buf:
                head, buf = buf.split(b'\r\n\r\n', 1)
                status = head.split(b'\r\n')[0].decode(errors='replace')
                if '200' not in status:
                    self.error = status
                    break
                chunked = b'chunked' in head.lower()
                header_done = True
            if header_done and buf:
                if chunked:
                    while True:
                        if b'\r\n' not in buf:
                            break
                        ln, rest = buf.split(b'\r\n', 1)
                        size = int(ln.split(b';')[0], 16) if ln else 0
                        if len(rest) < size + 2:
                            break
                        frames += rest[:size]
                        buf = rest[size + 2:]
                else:
                    frames += buf
                    buf = b''
                self.parse(frames)
            if not a.no_gga and time.time() - last_gga >= 5:
                try:
                    s.sendall(gga_sentence(lat, lon, h).encode())
                except OSError:
                    pass
                last_gga = time.time()
        s.close()

    def parse(self, frames):
        while len(frames) >= 6:
            if frames[0] != 0xD3:
                del frames[0]
                continue
            n = ((frames[1] & 3) << 8) | frames[2]
            if len(frames) < n + 6:
                return
            msg = bytes(frames[:n + 6])
            del frames[:n + 6]
            if crc24q(msg[:n + 3]) != int.from_bytes(msg[n + 3:n + 6], 'big'):
                self.crc_err += 1
                continue
            self.nbytes += len(msg)
            if self.a.out:
                self.raw += msg
            ty = getbits(msg, 24, 12)
            self.types[ty] = self.types.get(ty, 0) + 1
            sid = getbits(msg, 36, 12)
            self.staids.add(sid)
            if self.a.timeline:
                tw = getbits(msg, 48, 30) / 1000.0 if 1071 <= ty <= 1137 else -1.0
                self.timeline.append((time.time() - self.t_connect, ty, sid, tw))
            if 1071 <= ty <= 1137:
                tow = getbits(msg, 48, 30) / 1000.0
                if 1121 <= ty <= 1127:
                    tow += 14.0              # BDT -> GPST
                if self.first_msm is None:
                    self.first_msm = time.time() - self.t_connect
                age = (self.gps_now() - tow) % 604800.0
                if age > 302400:
                    age -= 604800
                sync = getbits(msg, 24 + 12 + 12 + 30, 1)
                if not sync:
                    self.ages.append(age)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--host', default='127.0.0.1')
    ap.add_argument('--port', type=int, default=2101)
    ap.add_argument('--mount', default='VRS_MSM7')
    ap.add_argument('--llh', type=float, nargs=3, required=True)
    ap.add_argument('--seconds', type=float, default=30)
    ap.add_argument('--clients', type=int, default=1)
    ap.add_argument('--spread', type=float, default=0.0, help='spread clients over +-km')
    ap.add_argument('--v2', action='store_true')
    ap.add_argument('--user'), ap.add_argument('--password')
    ap.add_argument('--no-gga', action='store_true')
    ap.add_argument('--out')
    ap.add_argument('--timeline', help='CSV of every message: client,t_s,type,staid,tow')
    ap.add_argument('--replay-start', type=float, help='GPS TOW of the replay start (replay mode)')
    ap.add_argument('--replay-offset', type=float, default=0.0)
    ap.add_argument('--speed', type=float, default=1.0)
    ap.add_argument('--replay-auto', action='store_true',
                    help='read the replay clock from the server /status (nrtk-svr rinex:// replay)')
    a = ap.parse_args()
    if a.replay_auto:
        st, t0 = http_status(a.host, a.port)
        if st.get('mode') == 'replay':
            a.replay_start, a.replay_offset = st['clock_tow'], 0.0
            a.replay_wall0 = t0
            print('replay clock: TOW %.1f' % a.replay_start)
    cl = [Client(a, i) for i in range(a.clients)]
    for c in cl:
        c.start()
        time.sleep(0.01)
    for c in cl:
        c.join()
    ages = [x for c in cl for x in c.ages]
    firsts = [c.first_msm for c in cl if c.first_msm is not None]
    errs = [c.error for c in cl if c.error]
    types = {}
    for c in cl:
        for k, v in c.types.items():
            types[k] = types.get(k, 0) + v
    print('clients %d, with data %d, errors %d %s' % (len(cl), len(firsts), len(errs), errs[:3]))
    print('message types:', dict(sorted(types.items())))
    print('station ids:', sorted(set().union(*[c.staids for c in cl]))[:10])
    print('CRC errors:', sum(c.crc_err for c in cl), ' bytes:', sum(c.nbytes for c in cl))
    if firsts:
        print('first MSM epoch after connect: mean %.2f s, max %.2f s' % (statistics.mean(firsts), max(firsts)))
    if ages:
        ages.sort()
        print('age of corrections: median %.2f s, p95 %.2f s, max %.2f s (%d epochs)' % (
            ages[len(ages) // 2], ages[int(0.95 * (len(ages) - 1))], ages[-1], len(ages)))
    if a.out and cl:
        open(a.out, 'wb').write(cl[0].raw)
        print('wrote', a.out)
    if a.timeline:
        with open(a.timeline, 'w') as f:
            f.write('client,t_s,type,staid,tow\n')
            for c in cl:
                for r in c.timeline:
                    f.write('%d,%.3f,%d,%d,%.3f\n' % ((c.idx,) + r))
        print('wrote', a.timeline)
    closed = sum(1 for c in cl if c.closed_by_server)
    if closed:
        print('connections closed by server:', closed)


if __name__ == '__main__':
    main()
