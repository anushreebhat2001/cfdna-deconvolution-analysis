import sys, time, subprocess, resource
t = time.time()
rc = subprocess.call(sys.argv[2:])
mem = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 / 1024
print(f"TIMING\t{sys.argv[1]}\t{time.time()-t:.1f}s\t{mem:.2f}GB\trc={rc}", flush=True)
