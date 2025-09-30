:loop
data > input.txt
brute < input.txt > answer.txt
test < input.txt > output.txt
fc output.txt answer.txt
if not errorlevel 1 goto loop
pause