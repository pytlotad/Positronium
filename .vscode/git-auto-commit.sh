#!/usr/bin/env bash
set -uo pipefail

repo_dir=$(git -C "$(dirname "$0")/.." rev-parse --show-toplevel) || exit 1
cd "$repo_dir" || exit 1

if ! command -v inotifywait >/dev/null 2>&1; then
    echo "ERROR: inotifywait is required. Install inotify-tools." >&2
    exit 1
fi

commit_changes() {
    # Editors often save several related files a moment apart. Group them into
    # one coherent commit instead of producing one commit per filesystem event.
    sleep 1
    git add -A
    git diff --cached --quiet && return 0

    local changed_files summary body
    changed_files=$(git diff --cached --name-only)
    summary="Aktualizuj projekt"
    body=""

    if grep -qx 'positronium.cpp' <<<"$changed_files"; then
        summary="Aktualizuj symulację pozytonium"
        body+="- Zmień implementację i zachowanie symulacji."$'\n'
    fi
    if grep -qx 'README.md' <<<"$changed_files"; then
        body+="- Zaktualizuj dokumentację projektu."$'\n'
    fi
    if grep -qx 'Makefile' <<<"$changed_files"; then
        body+="- Dostosuj budowanie i uruchamianie programu."$'\n'
    fi
    if grep -q '^\.vscode/' <<<"$changed_files"; then
        if [[ "$summary" == "Aktualizuj projekt" ]]; then
            summary="Skonfiguruj automatyczne commity w VS Code"
        fi
        body+="- Zaktualizuj konfigurację środowiska VS Code."$'\n'
    fi
    if grep -qx '.gitignore' <<<"$changed_files"; then
        body+="- Zaktualizuj reguły ignorowania plików roboczych."$'\n'
    fi

    body+=$'\nPliki:\n'
    while IFS= read -r changed_file; do
        body+="- ${changed_file}"$'\n'
    done <<<"$changed_files"

    if git commit -m "$summary" -m "$body"; then
        echo "Committed: $summary"
        sync_changes
    else
        echo "WARNING: automatic commit failed; changes remain staged." >&2
    fi
}

sync_changes() {
    local branch
    branch=$(git branch --show-current)
    if [[ -z "$branch" ]]; then
        echo "WARNING: automatic sync skipped in detached HEAD state." >&2
        return 0
    fi
    if ! git remote get-url origin >/dev/null 2>&1; then
        echo "WARNING: automatic sync skipped: remote 'origin' is not configured." >&2
        return 0
    fi

    echo "Syncing branch: $branch"
    if git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}' >/dev/null 2>&1; then
        if ! GIT_TERMINAL_PROMPT=0 git pull --rebase --autostash; then
            echo "WARNING: automatic pull failed; local commit is preserved." >&2
            return 0
        fi
        if ! GIT_TERMINAL_PROMPT=0 git push; then
            echo "WARNING: automatic push failed; local commit is preserved." >&2
            return 0
        fi
    elif ! GIT_TERMINAL_PROMPT=0 git push --set-upstream origin "$branch"; then
        echo "WARNING: automatic push failed; local commit is preserved." >&2
        return 0
    fi
    echo "Sync complete: $branch"
}

echo "Watching repository for saved files..."
inotifywait -m -r -e close_write,move,create,delete \
    --exclude '(^|/)(\.git|\.positronium-frames\.cache|\..*~?|.*\.swp$|.*\.tmp$)' \
    --format '%w%f' . 2>&1 | while IFS= read -r changed_path; do
        case "$changed_path" in
            "Setting up watches."*|"Watches established.")
                echo "$changed_path"
                ;;
            *)
                echo "Detected change: $changed_path"
                commit_changes
                ;;
        esac
    done
