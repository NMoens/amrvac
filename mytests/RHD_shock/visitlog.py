# Visit 2.13.2 log file
ScriptVersion = "2.13.2"
if ScriptVersion != Version():
    print "This script is for VisIt %s. It may not work with version %s" % (ScriptVersion, Version())
ShowAllWindows()
Source("/lhome/nicolasm/Visit/2.13.2/linux-x86_64/bin/makemovie.py")
OpenComputeEngine("altair.ster.kuleuven.be", "-forcestatic")
RestoreSession("/lhome/nicolasm/amrvac/mytests/RHD_shock/movie_Er.session", 0)
