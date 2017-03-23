import sys




def readData(data,fname):
    j=0;
    with open(fname) as f:
        content = f.readlines()
        for line in content:
            elems = line.split()
            for elem in elems:
                data[j] = float(elem)
                j = j+1



def calcDiff(data1,data2):
    print "[Element] : [diff value] - [diff relative %]  "
    val = ""
    for i in xrange(1,len(data1)):
        val = str(i) + ": " + str(abs(data1[i] - data2[i])) + " "
        calc = 0 
        if (data1[i]!=0):
            calc = 100*abs(data1[i] - data2[i])/data1[i]
            val = val +  str(abs(data1[i] - data2[i])/data1[i]) + " "
            if (calc > 5):
                val = val + " <<<<<<===== WARNING over 5% diff: " + str(abs(data1[i] - data2[i]))
            print  val
        val = ""
        if (i%74==0):
            print  "-------------------"


data1 = [0 for y in range(400)] 
data2 = [0 for y in range(400)] 

def calcError(prefix1,prefix2):

    readData(data1,prefix1)
    readData(data2,prefix2)
    calcDiff(data1,data2)



## MAIN

fname1 = sys.argv[1]
fname2 = sys.argv[2]

calcError(fname1,fname2)




