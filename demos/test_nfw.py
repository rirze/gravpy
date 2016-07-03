from gravpy.models import nfwf

print "running nfwf module"
e = 0.0001
print nfwf.single_eval(0.100001,0.0999987,(1,0,0,e,0,0.1))[1:]
print nfwf.single_eval(0.1,0.1,(1,0,0,0,0,0.1))[1:]

