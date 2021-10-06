import os

print("""
	From or to local machine?
	""")
whichWay = str(input ("Type 'from' or 'to': "))
print("""
	Are you copying a directory or file?
	""")
isDir = str(input ("Type 'dir' or 'file': "))

arg = ' -r ' if isDir == 'dir' else ' '
remote = 'jmccord9@login-phoenix.pace.gatech.edu:~/p-amedford6-0/rich_project_chbe-medford/p-dopants/fall2021/'
local = '~/Dropbox\ \(GaTech\)/pBlock/fall2021/'
scp = 'scp' + arg

print("\n        Input local directory or file after "+local+"\n")
whichLocal = str(input ("Local directory or file: " ))
print("\n        Input remote directory or file after "+remote+"\n")
whichRemote = str(input ("Remote directory or file: " ))
print("""
	I will execute the following command. Are you sure you would like to continue?
	""")

local += whichLocal
remote += whichRemote

if whichWay == 'from':
	command = scp+local+" "+remote
elif whichWay == 'to':
	command = scp+remote+" "+local

print("\n"+command+"\n")
procede = str(input ("Type 'yes' or 'no': "))
print("Executing command!") if procede == 'yes' else None
os.system(command) if procede == 'yes' else print('Exiting!')
# os.system("echo This is a temp command.") if procede == 'yes' else print('Exiting!')
