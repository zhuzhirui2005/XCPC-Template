# chmod +x duipai.sh

SCRIPT="$(readlink -f "$0")"
if [ -z "$TERM" ] || [ "$TERM" = "dumb" ]; then
  if command -v x-terminal-emulator >/dev/null 2>&1; then
    exec x-terminal-emulator -e bash -lc "\"$SCRIPT\""
  elif command -v gnome-terminal >/dev/null 2>&1; then
    exec gnome-terminal -- bash -lc "\"$SCRIPT\""
  elif command -v konsole >/dev/null 2>&1; then
    exec konsole -e bash -lc "\"$SCRIPT\""
  fi
fi
cd -- "$(dirname -- "$SCRIPT")" || exit 1
while true; do
  ./data > input.txt || break
  ./brute < input.txt > answer.txt || break
  ./test < input.txt > output.txt || break
  diff -q output.txt answer.txt >/dev/null 2>&1 && continue || break
done
read -n 1 -s -r -p "Press any key to continue..."
echo