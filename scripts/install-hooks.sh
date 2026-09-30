#!/usr/bin/env bash
# Opt-in: installs a pre-push hook that runs scripts/check.sh.
set -euo pipefail
cd "$(dirname "$0")/.."
hook=.git/hooks/pre-push
printf '#!/usr/bin/env bash\nexec "$(git rev-parse --show-toplevel)/scripts/check.sh"\n' > "$hook"
chmod +x "$hook"
echo "installed $hook"
