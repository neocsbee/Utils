"""
Created: Jan 20, 2013: (use this file only after unzipping the dropbox file)
start from the directory that contains all the studentfolders from dropbox
from lms. This folder should not have any other directories other than
student submissions, but other files are allowed. 

Modified: Jan 21, 2012: finding out multiple submissions from students, based on
student last and first name

Result of running this script: messages from students displayed to stdout
(might want to redirect them to a file to read offline) and contents of zip
are taken and added to the student directory without creating any additional directory
structure. This script is capable of going deep inside the directory structure inside the zip file
and just copies the files along the way to the parent directory

Jan 23, 2013: did some modification to save the .m files with the student
name attached to them to the parent directory that way all the files are in
same place and directory switch could be avoided. May not be able to fill all
rows corresponding to a student at once.

Feb 2, 2013: Further improvements to this file, to write the names alone of .m
files to a m file for batch execution of scripts. Since there are two target
files to be written, had to close the source and reopen it to set the file
pointer at offset 0 inside the file.
    
Got the code from :
    stackoverflow.com/questions/4917284/extract-files-from-zip-without-keeping-the-structure-using-python-zipfile.
    The second answer also looks good and is without shutil usage and
    straightforward to understand
For file extension check:
    stackoverflow.com/questions/5899497/checking-file-extension

Feb 15, 2013: Major changes: 
    1) Similar to hw2, write the script files with student names attached to the parent
    directory. As file name overrides function names, it is possible to do this. 
    2) Cannot write function calls to a mscript because I need to pass arguments to 
    functions and they are different for each function.

Feb 22, 2013: Now the mscript is written for executing scripts/functions in alphabetical
order, helps to reduce the up/down movement in gradeentry.ods file. Logic to identify
a file to be function or script is finalized. This homework in particular has .doc based
additional file submission.
"""
import os
import shutil
import zipfile
import subprocess
import re
parentdir = './'
sodirs = set()
mfilecount = 0
ofilecount = 0
#rename the student folders such that last name is listed first and
#first name is listed second
dirlist = os.listdir(parentdir)
for x in dirlist:
    if x.find(' ')!=-1:
        xtemp = x.split(' ')
        if len(xtemp) > 2:
            y = xtemp[1]+' '+xtemp[0]+' '+' '.join(xtemp[2:])
        else:
            y = xtemp[1]+' '+xtemp[0]
        os.rename(x, y)
#sort the names so that entering grades in gradeentry.ods is easier as the names in
#gradeentry.ods is sorted according to last names and so is grade.htm file
dirlist = sorted(os.listdir(parentdir))
fm = open('mscript.m','w');
for x in dirlist:
    if os.path.isdir(x):
        #extract the first two strings, join them as one, check their presence in set
        #of already seen directories - this checks if a student name is redundant
        key = ' '.join(x.split()[:2])
        if key in sodirs:
            print 'multiple submission:', x
        else:
            sodirs.add(key)
        for y in os.listdir(parentdir+x):
            ext = os.path.splitext(y)[-1].lower()
            if  ext not in ('.zip','.htm'):
                    print 'see', x, 'wrong ext', ext
                    continue
            if ext == '.htm':
                htmfile = parentdir+x+'/'+y
                statinfo = os.stat(htmfile)
                if statinfo.st_size!=0:
                    print 'Message from '+x
                    print subprocess.check_output(["cat", htmfile])
            if ext == '.zip':
                print>>fm, "cd('{0}')".format(parentdir+x)
                archfile = parentdir+x+'/'+y
                #print 'archfile:', archfile
                zf=zipfile.ZipFile(archfile,'r')
                for member in zf.namelist():
                    filename=os.path.basename(member)
                    mname = os.path.splitext(filename)[0]
                    mext = os.path.splitext(filename)[-1]
                    if not filename:
                            continue
                    if mext.lower() != '.m' and mext.lower()!='.m~' and \
                            mext.lower()!='.doc':#handle docx too
                            ofilecount = ofilecount + 1
                            continue
                    try:
                        if mext == '.m~':
                            filename = mname+".m"
                            mext = ".m"
                        if mext=='.m':
                            mfilecount = mfilecount + 1
                            flag =  re.match(r'^[a-zA-Z][\w]*',mname)
                            if flag == None: #there is no match
                                print "{0} : wrong file name {1}"\
                                        .format(x, filename)
                        source = zf.open(member)
                        stuname = ''.join(x.split()[:2])
                        #wb flag : write binary,overwrite if already exists
                        #first target is inside the student directory
                        target1 = file(os.path.join(parentdir+x, filename), "wb")
                        shutil.copyfileobj(source,target1)
                        source.close()
                        target1.close()
                        #based on the file content write as script or function
                        funcfound = False
                        if mext==".m":
                            filecont = subprocess.check_output(["cat",\
                                    os.path.join(parentdir+x,filename)])
                            if filecont.count("function")!=0:
                                mylist = filecont.splitlines()
                                for line in mylist:
                                    funcidx = line.find('function')
                                    eqidx = line.find('=')
                                    opidx = line.find('(')
                                    epidx = line.find(')')
                                    if funcidx != -1 and \
                                            eqidx != -1 and opidx != -1 \
                                            and epidx !=-1:
                                                if funcidx < eqidx and \
                                                        eqidx < opidx and\
                                                        opidx < epidx:
                                                    funcfound = True
                                                    print>>fm,"{0}( )".format(mname)
                                                    print>>fm,"{0}( )".format(mname)
                                                    print>>fm,"{0}( )".format(mname)
                                                    print>>fm,"clear\nclc"
                                                    print>>fm, \
                                                    "input('grade {0} {1}');"\
                                                    .format(stuname, mname)
                                                    break
                            if not funcfound:
                                print>>fm,"{0}".format(mname)
                                print>>fm, "input('grade {0} {1}');".format(stuname, mname)
                                print>>fm, "clear\nclc"
                                print>>fm,"{0}".format(mname)
                                print>>fm, "input('grade {0} {1}');".format(stuname, mname)
                                print>>fm, "clear\nclc"
                                print>>fm,"{0}".format(mname)
                                print>>fm, "input('grade {0} {1}');".format(stuname, mname)
                                print>>fm, "clear\nclc"
                    except zipfile.BadZipfile:
                        print 'bad zip file in the folder:', x
                        print 'skipping the copy process'
                zf.close()
        #one student directory grading is over, come out to parent directory
        print>>fm, "cd .."
        print '{0} : m-files: {1}, other file count:{2}'\
                .format(x, mfilecount, ofilecount)
        mfilecount = 0
        ofilecount = 0
fm.close()
del sodirs
