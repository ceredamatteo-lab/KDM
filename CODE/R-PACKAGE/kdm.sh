#!/usr/bin/env bash
# kdm.sh - start the kdm R container with host X11 and extra mounts
# (the current directory is always mounted at its own path and used as start dir)
#
# Usage:
#   ./kdm.sh [mount ...] [-- command ...]
#
# Each mount is a host path, mounted at the same absolute path in the
# container. Append :ro for read-only (e.g. /data/refs:ro).
#
# Examples:
#   ./kdm.sh
#   ./kdm.sh /data/genomes ~/refs:ro
#   ./kdm.sh /data/genomes -- Rscript analysis.R
#
# Environment overrides: KDM_IMAGE (default kdm),
#                        KDM_GROUP (default: your primary group on the host)

set -euo pipefail

IMAGE="${KDM_IMAGE:-kdm}"
GROUP_NAME="${KDM_GROUP:-}"

die() { echo "kdm.sh: $*" >&2; exit 1; }

# --- user / group ----------------------------------------------------------
uid=$(id -u)
if [[ -n "$GROUP_NAME" ]]; then
    gid=$(getent group "$GROUP_NAME" | cut -d: -f3 || true)
    [[ -n "$gid" ]] || die "group '$GROUP_NAME' (KDM_GROUP) not found"
else
    gid=$(id -g)
fi

# -t only when attached to a terminal, so the script also works in pipes/cron
if [[ -t 0 && -t 1 ]]; then tty_flags=(-it); else tty_flags=(-i); fi

# current directory: mounted at its own (resolved) path and used as start dir
cwd=$(pwd -P)

args=( "${tty_flags[@]}" --rm
       --user "$uid:$gid"
       -e HOME=/tmp
       -v "$cwd:$cwd"
       -w "$cwd" )

# all the host user's supplementary groups (numeric, so they need not
# exist inside the image); the primary group above decides new files' group
for g in $(id -G); do
    [[ "$g" == "$gid" ]] || args+=( --group-add "$g" )
done

# --- X11 (host display) ----------------------------------------------------
if [[ -n "${DISPLAY:-}" ]]; then
    args+=( -e "DISPLAY=$DISPLAY" -v /tmp/.X11-unix:/tmp/.X11-unix )
    xauth_file="${XAUTHORITY:-$HOME/.Xauthority}"
    if [[ -f "$xauth_file" ]]; then
        args+=( -v "$xauth_file:/tmp/.Xauthority:ro" -e XAUTHORITY=/tmp/.Xauthority )
    fi
    # ssh -X gives a TCP display such as localhost:10.0 -> needs host network
    if [[ "$DISPLAY" != :* ]]; then
        args+=( --net=host )
    fi
else
    echo "kdm.sh: DISPLAY is not set, X11 plot windows will not be available" >&2
fi

# --- extra mounts ----------------------------------------------------------
while (( $# )); do
    case "$1" in
        --) shift; break ;;
        -h|--help) sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*) die "unknown option: $1 (use -- before the command to run)" ;;
        *)
            src=$1; shift
            mode=""
            if [[ "$src" == *:ro ]]; then
                mode=":ro"
                src=${src%:ro}
            fi
            [[ "$src" != *:* ]] || die "paths containing ':' are not supported: $src"
            [[ -e "$src" ]] || die "mount source does not exist: $src"
            src=$(realpath "$src")
            # the current directory is already mounted; Docker rejects duplicates
            [[ "$src" == "$cwd" ]] && continue

            args+=( -v "$src:$src$mode" )
            ;;
    esac
done

# --- user and group names inside the container -----------------------------
# Start from the image's own /etc/passwd and /etc/group and append the host
# user and their groups (looked up with getent, so LDAP/SSSD groups work too).
# The result is bind-mounted read-only over the image's files.
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

docker run --rm --entrypoint sh "$IMAGE" -c 'cat /etc/passwd; echo ===; cat /etc/group' \
  | awk -v d="$tmpdir" '/^===$/ {f = d "/group"; next} {print > (f ? f : d "/passwd")}'

has_id() { awk -F: -v id="$2" '$3 == id {f = 1} END {exit !f}' "$1"; }
drop_id() { awk -F: -v id="$2" '$3 != id' "$1" > "$1.new" && mv "$1.new" "$1"; }

user=$(id -un)
for g in "$gid" $(id -G); do
    grep -q ":x:$g:$user\$" "$tmpdir/group" && continue   # already added
    # ids >= 1000 belong to ordinary accounts in the image (Ubuntu 24.04 ships
    # 'ubuntu' as 1000): the host name wins. System groups below 1000 are kept.
    if has_id "$tmpdir/group" "$g"; then
        (( g >= 1000 )) || continue
        drop_id "$tmpdir/group" "$g"
    fi
    name=$(getent group "$g" | cut -d: -f1 || true)
    [[ -n "$name" ]] && echo "$name:x:$g:$user" >> "$tmpdir/group"
done
drop_id "$tmpdir/passwd" "$uid"
echo "$user:x:$uid:$gid::/tmp:/bin/bash" >> "$tmpdir/passwd"

args+=( -v "$tmpdir/passwd:/etc/passwd:ro" -v "$tmpdir/group:/etc/group:ro" )

# not exec: the script must stay alive to remove the temporary files afterwards
docker run "${args[@]}" "$IMAGE" "$@"
