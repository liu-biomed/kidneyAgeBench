"""Parallel HTTPS range downloader for NCBI GEO files.
The local route throttles each HTTPS connection (~4 KB/s) and corrupts FTP streams, so each file
is fetched as many byte ranges over concurrent connections, reassembled, and integrity-checked."""
import sys, os, gzip, zipfile, tarfile, threading, queue, time, urllib.request
opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
CH = 192 * 1024
def size_of(url):
    r = opener.open(urllib.request.Request(url, method="HEAD", headers={"User-Agent": "curl/8"}), timeout=60)
    return int(r.headers["Content-Length"])
def fetch(url, a, b, path):
    have = os.path.getsize(path) if os.path.exists(path) else 0
    need = b - a + 1
    while have < need:
        try:
            req = urllib.request.Request(url, headers={"Range": f"bytes={a+have}-{b}", "User-Agent": "curl/8"})
            with opener.open(req, timeout=120) as r, open(path, "ab") as f:
                while True:
                    buf = r.read(16384)
                    if not buf: break
                    f.write(buf)
        except Exception:
            time.sleep(2)
        have = os.path.getsize(path) if os.path.exists(path) else 0
    if have > need:  # defensive: truncate overshoot
        with open(path, "r+b") as f: f.truncate(need)
def check(path):
    try:
        if path.endswith(".gz"):
            with gzip.open(path, "rb") as f:
                while f.read(1 << 20): pass
        elif path.endswith(".xlsx"): zipfile.ZipFile(path).testzip()
        elif path.endswith(".tar"): tarfile.open(path).getmembers()
        return True
    except Exception as e:
        print("CHECK FAIL", path, e, flush=True); return False
def main(listfile, nthreads=64):
    jobs = [l.split() for l in open(listfile) if l.strip()]
    q = queue.Queue(); plan = []
    for url, out in jobs:
        if os.path.exists(out) and check(out): print("ok", out, flush=True); continue
        if os.path.exists(out): os.remove(out)
        S = size_of(url); parts = []
        for i, a in enumerate(range(0, S, CH)):
            p = f"{out}.part{i:05d}"; parts.append(p); q.put((url, a, min(a + CH, S) - 1, p))
        plan.append((out, parts, S)); print("queued", out, S, len(parts), flush=True)
    def worker():
        while True:
            try: job = q.get_nowait()
            except queue.Empty: return
            fetch(*job); q.task_done()
    ts = [threading.Thread(target=worker, daemon=True) for _ in range(nthreads)]
    [t.start() for t in ts]; [t.join() for t in ts]
    for out, parts, S in plan:
        with open(out, "wb") as f:
            for p in parts:
                with open(p, "rb") as g: f.write(g.read())
                os.remove(p)
        print("done", out, os.path.getsize(out) == S, check(out), flush=True)
if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 64)
