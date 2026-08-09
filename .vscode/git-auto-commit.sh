#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.." || exit 1
if ! command -v inotifywait >/dev/null 2>&1; then
  echo "ERROR: inotifywait is required. Install inotify-tools." >&2
  exit 1
fi

echo "Watching repository for saved files..."

while inotifywait -m -r -e close_write --exclude '(^|/)(\.git|\.vscode|\..*~?|.*\.swp$|.*\.tmp$)' . | while read -r path action file; do
    echo "Detected save: ${path}${file}"
    git add -A
    if git diff --cached --quiet; then
        continue
    fi
    files=$(git diff --cached --name-only | tr '\n' ' ' | sed 's/ $//')
    summary="Auto commit on save"
    if git diff --cached --name-only | grep -q '^\.vscode/git-auto-commit\.sh$'; then
        summary="Add auto-commit watcher"
    fi
    if git diff --cached --name-only | grep -q '^\.vscode/tasks\.json$'; then
        summary="Add VS Code auto-commit task"
    fi
    if git diff --cached --name-only | grep -q '^\.vscode/settings\.json$'; then
        summary="Update editor auto-save / git settings"
    fi
    if git diff --cached --name-only | grep -q '^positronium\.cpp$'; then
        diff=$(git diff --cached --unified=0 -- positronium.cpp)
        if echo "$diff" | grep -q 'SetLineColor(kAzure' && echo "$diff" | grep -q 'SetLineColor(kRed'; then
            summary="Color electron blue and positron red"
        elif echo "$diff" | grep -q 'labelFor(' && echo "$diff" | grep -q 'deltaLabelFor('; then
            summary="Update energy readouts"
        elif echo "$diff" | grep -q 'git-auto-commit' ; then
            summary="Update auto-commit integration"
        else
            summary="Update positronium code"
        fi
    fi
    timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    commit_message="${summary} (${timestamp}): ${files}"
    git commit -m "$commit_message"
    echo "Committed: $commit_message"
done
done
