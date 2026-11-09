#!/bin/bash
# diagnostics: show command lines / fds of running qemu-system-arm processes, optionally kill them (arg: kill)
for p in $(pgrep -x qemu-system-arm); do
    echo "== pid $p  (etime $(ps -o etime= -p $p)) =="
    tr '\0' ' ' < /proc/$p/cmdline; echo
    ls -la /proc/$p/fd 2>/dev/null | awk '{print $9, $10, $11}' | head -6
    if [ "$1" = "kill" ]; then kill -9 $p; echo "killed $p"; fi
done
for p in $(pgrep -x timeout); do [ "$1" = "kill" ] && kill -9 $p 2>/dev/null; done
exit 0
