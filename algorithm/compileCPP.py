#! /usr/bin/env python

import os

a = raw_input('Pfad: ')
b = raw_input('Name: ')

if a:
    path = a + '/' + b
else:
    path = b

version = raw_input('version')
if version == '2.7':
    cmd1 = 'g++ -c -fPIC -I/usr/include/python2.7 ' + path + '.cpp'
    cmd2 = 'g++ -shared ' + path + '.o -o '+ path + '.so -lpython2.7 -lboost_python'
elif version == '3.3':
    print "yippie"
    cmd1 = 'g++ -c -fPIC -I/usr/include/python3.3 ' + path + '.cpp'
    cmd2 = 'g++ -shared ' + path + '.o -o '+path + '.so -lpython3.3 -lboost_python'
else:
    cmd1 = 'g++ -c -fPIC -I/usr/include/python2.6 ' + path + '.cpp'
    cmd2 = 'g++ -shared ' + path + '.o -o '+ path + '.so -lpython2.6 -lboost_python'

print cmd1
os.system(cmd1)
print cmd2
os.system(cmd2)
raw_input("Press enter to quit.")

