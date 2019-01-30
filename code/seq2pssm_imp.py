import sys, getopt
import numpy
import io
import os
import shlex
import subprocess
import numpy
import glob
currdir = os.getcwd();
def main(argv):
   inputfile = ''
   outputfile = ''
   database = ''
   
   if len(argv[1:]) == 0:
#        print ('seq2motif.py -i <inputfile> -o <outputfile> -w <window> -x <extension>')
        print ("\nProgram for computing PSSM matrix in column format without any normalization.\n")
        print ("Usage:python seq2pssm_imp.py -i inputFile -o outputFile -d database\n")
        print ('-i\tInputFile (SFASTA FORMAT)\n-o\toutputFile\n-d\tDatabase\n')
        sys.exit()
   try:
       # opts is a list of returning key-value pairs, args is the options left after striped
       # the short options 'hi:o:', if an option requires an input, it should be followed by a ":"
       # the long options 'ifile=' is an option that requires an input, followed by a "="
      opts, args = getopt.getopt(argv,"i:o:d:",["ifile=","ofile=","dbase="])
   except getopt.GetoptError:
      print ('Error')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('seq2pssm_imp.py -i <inputfile> -o <outputfile> -d <database>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      elif opt in ("-d", "--dbase"):
         database = sys.argv[6]	 
   #S = currdir+'/blastpr/blastpgp -d '+currdir+'blastpr/data/'+database+' -i '+inputfile+' -j 3 -C ',outputfile;
#S = '/home/users/piyush/pssm_run/blastpr/blastpgp -d /home/gpsr/data/blastdata/{} -i {} -j {} -C {}'.format(database, inputfile, 3, outputfile)
   S = '{}/blastpr/blastpgp -d {}/blastpr/data/{} -i {} -j {} -C {}'.format(currdir, currdir, database, inputfile, 3, outputfile)
   print('###########\t',S);
   os.system(S)

   filename, file_extension = os.path.splitext(inputfile)
   filename1, file_extension_o = os.path.splitext(outputfile)
   filename_o, file_extension_o = os.path.splitext(outputfile) 

   os.rename(outputfile, filename_o + '.chk')
   outputfile1 = filename_o + ".chk"
   temp1 =filename + ".sn"
   #inputfile > temp1
     
   C2 = 'echo {} > {}'.format(inputfile,temp1)
   os.system(C2)

   temp2 = filename_o + ".pn"
   #outputfile > temp2     
   C1 = 'echo {} > {}'.format(outputfile1,temp2)
   os.system(C1)
  # P = currdir+'/blastpr/makemat -P '+filename;
   P = '{}/blastpr/makemat -P {}'.format(currdir,filename)
   #P = '/home/users/piyush/pssm_run/blastpr/makemat -P {}'.format(filename)
   print('********',P);
   os.system(P)		
   
   dir = '.'
   allfiles = os.listdir(dir)
   files = [ fname for fname in allfiles if fname.endswith('.mtx')]
   print(files)
   x = 0
   while x < len(files) :
    if files[x] == filename + ".mtx" :
     mtx_file = files[x]
    x += 1

   orig_stdout = sys.stdout
   n1 = open('PSSMProfile.txt','w')     
   sys.stdout = n1
   
   numpy.set_printoptions(threshold='nan')
   with open(mtx_file,'r') as f:
       g = list(f)   
	
   count = 0
   i = 0
   j = 0
   i1 = 0
   j1 = 0
   c1 = 0
   c2 = 0
   q1 = 0
   q2 = 0
   l1 =0
   s1 = 0  
   r = 0
   ct1 = 0   
   x1 = 0
   y1 = 0
   z1 = 0
   x2 = 0
   y2 = 0
   z2 = 0
   l = 0
   w = 0
   d1 = 0
   d2 = 0
   d3 = 0
   while i1 < len(g) :
    g[i1] = g[i1].replace('\n','')
    i1 = i1 + 1 
	
   while i < len(g):
    if 'e' in g[i] :
     j = i
    else :
     count += 1
    i += 1	
	
   k = j - (len(g)-count)+1
   g[k:j+1] = []	
   p1 = k
   
   while p1 < len(g) :
    g[p1] = g[p1].split(" ")
    p1 += 1
	
   while s1 < len(g[2]) :
    if g[2][s1] == '' :
     c1 += 1
    s1 += 1	
	
   hh =  len(g[2])-c1

   bb = numpy.zeros(shape=(len(g)-2,hh))	
   while q2 < len(g) :
    if q2 > k -1 :
     while q1 < len(g[2]) :
      if g[q2][q1] == '':
       c1 += 1
      else :
       bb[q2-k][c2] = g[q2][q1]
       c2 += 1
      q1 += 1
     q1 = 0
    c2 = 0
    q2 += 1
	
   while r < len(bb[0]) :
    if bb[0][r] == -32768. :
     ct1 = ct1 + 1
    r += 1	
	
   ct2 = len(bb[0]) - ct1 
   cc = numpy.zeros(shape=(len(bb),ct2))
   dd = numpy.zeros(shape=(len(bb),ct1))	
   
   while x2 < len(bb) :
    while y2 < len(bb[0]) :
     if bb[x2][y2] != -32768. :
      cc[x2][z2] = bb[x2][y2]
      z2 += 1
     y2 += 1
    z2 = 0
    y2 = 0
    x2 += 1
	
   ee = numpy.zeros(shape=(len(bb),ct2+1), dtype = object)
   while d1 < len(ee) :
    while d2 < ct2 + 1 :
     ee[d1][0] = g[1][d1] 
     if d2 > 0 :
      ee[d1][d2] = cc[d1][d2 - 1]  
     d2 += 1
    d2 = 0
    d1 += 1 

   ff = ee.tolist()
   print(str(ff).replace('[','').replace('], ','\n').replace(']]',''))
   n1.truncate()

main(sys.argv[1:])   
   
