#!/usr/bin/env bash
# Stage exactly the paths listed in git-files.txt, and nothing else.
#
#   ./git-add-listed.sh --dry     show what would be staged, change nothing
#   ./git-add-listed.sh           stage them
#
# Why this exists rather than a bare `git add --pathspec-from-file=git-files.txt`:
# that form rejects comment and blank lines with "empty string is not a valid
# pathspec", so the manifest could not be annotated. Here the manifest is filtered
# and handed to git on stdin instead.
#
# The manifest SELECTS; .gitignore still VETOES. A listed path that .gitignore
# excludes is reported and skipped rather than added, which is what keeps
# zotero_key.txt and the data trees out even if a line for one is added by
# mistake. Nothing here uses `git add -f`, deliberately.
set -euo pipefail

cd "$(dirname "$0")"
MANIFEST=git-files.txt

command -v git >/dev/null || { echo "git not on PATH; run: ml load system git/2.45.1" >&2; exit 1; }
case "$(git --version)" in
  *" 1."*) echo "git $(git --version) is too old for --pathspec-from-file." >&2
           echo "Run: ml load system git/2.45.1" >&2; exit 1;;
esac
[ -f "$MANIFEST" ] || { echo "$MANIFEST not found" >&2; exit 1; }

# strip comments and blank lines
paths=$(grep -vE '^[[:space:]]*(#|$)' "$MANIFEST" || true)
[ -n "$paths" ] || { echo "$MANIFEST lists no paths" >&2; exit 1; }
n=$(printf '%s\n' "$paths" | wc -l)

if [ "${1:-}" = "--dry" ]; then
    echo "would stage $n path(s) from $MANIFEST:"
    printf '%s\n' "$paths" | git add --pathspec-from-file=- --dry-run
    exit 0
fi

printf '%s\n' "$paths" | git add --pathspec-from-file=-
echo "staged from $MANIFEST ($n path(s) listed)"
echo "now in the index: $(git diff --cached --name-only | wc -l) file(s)"
echo
echo "review with : git status --short"
echo "commit with : git commit -m 'your message'"
