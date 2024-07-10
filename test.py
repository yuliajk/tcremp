import os
os.remove('test.txt') if os.path.exists('test.txt') else print('no file')