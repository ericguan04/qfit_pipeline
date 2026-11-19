# Stage state from the marker files. Sourced, not executed.

STAGES="pre_qfit qfit post_qfit"

marker_value() {  # marker_value FILE KEY
  sed -n "s/^$2=//p" "$1" 2>/dev/null | head -n 1
}

# True if the job behind a .running marker still exists. A SLURM job is
# looked up with squeue; a local one by pid, if it ran on this host.
# Unknown (another host, no squeue) counts as alive.
running_alive() {
  local job host pid
  job=$(marker_value "$1" job)
  host=$(marker_value "$1" host)
  pid=$(marker_value "$1" pid)
  if [ -n "${job}" ] && command -v squeue >/dev/null 2>&1; then
    [ -n "$(squeue -h -j "${job}" 2>/dev/null)" ]
  elif [ -z "${job}" ] && [ "${host}" = "$(hostname)" ] && [ -n "${pid}" ]; then
    kill -0 "${pid}" 2>/dev/null
  else
    return 0
  fi
}

# stage_state STRUCTURE_DIR STAGE -> done | failed | running | interrupted | pending
stage_state() {
  local s="$1/status/$2"
  if [ -f "${s}.done" ]; then echo done
  elif [ -f "${s}.running" ]; then
    if running_alive "${s}.running"; then echo running; else echo interrupted; fi
  elif [ -f "${s}.failed" ]; then echo failed
  else echo pending
  fi
}
