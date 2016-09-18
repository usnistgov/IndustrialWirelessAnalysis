
import os

urlbase = 'https://s3.amazonaws.com/nist-ics-wireless/measurements/'

myDirPaths = [
	"active_meas/GBurg", 
	"active_meas/AAplant", 
	"active_meas/Boulder", 
	"active_meas/Boulder_c",
	"passive_meas/NIST MachineShop", 
	"passive_meas/TransmissionAssembly/TX1 position M18", 
	"passive_meas/TransmissionAssembly/TX2 position P10_Q10"]
print(myDirPaths)
myOutFiles = [
	"gburg_active.htm", 
	"aaplant_active.htm", 
	"boulder_active.htm", 
	"boulderc_active.htm",
	"gburg_passive.htm", 
	"aaplant_passive1.htm", 
	"aaplant_passive2.htm" ]
print(myOutFiles)

'''
# This section lists groups into separate listings
ndirs = len(myDirPaths)
for ii in range(ndirs):
	print(ii)
	myDirPath = myDirPaths[ii]
	myOutFile = myOutFiles[ii]
	print("Working on %s ==> %s" % (myDirPath, myOutFile))
	myFiles = os.listdir(myDirPath)
	f = open(myOutFile, "w")
	f.write('<!DOCTYPE htm PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">')
	f.write('<htm>\n<head>\n</head>\n<meta http-equiv="Content-Type" content="text/htm; charset=utf-8">\n<body>\n<ul>\n')
	f.writelines(['<li><a href="%s%s/%s">%s</a></li>\n' % (urlbase, myDirPath, f, f) for f in myFiles if f.endswith('.mat') or f.endswith('.csv') or f.endswith('.txt') or f.endswith('.tar')])
	f.write('</ul>\n</body>\n</htm>\n')
	f.close()
'''

# This section lists groups into a single file
ndirs = len(myDirPaths)
myOutFile = "downloads.htm"
f = open(myOutFile, "w")
myGroupHeaders = [
	"Gaithersburg Machine Shop ACTIVE Mobile Measurements", 
	"Automotive ACTIVE Mobile Measurements", 
	"Boulder Steam Plant ACTIVE Mobile Measurements", 
	"Boulder Steam Plant ACTIVE Stationary Measurements",
	"Gaithersburg Machine Shop PASSIVE Measurements", 
	"Automotive PASSIVE Measurements TX1 Location", 
	"Automotive PASSIVE Measurements TX2 Location" ]
anchors = [
	"gburgm",
	"autom",
	"boulderm",
	"boulderc",
	"gburgp",
	"autop1",
	"autop2"
]
f.write('<!DOCTYPE htm PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">')
f.write('<htm>\n<head>\n')
  
f.write('</head>\n<meta http-equiv="Content-Type" content="text/htm; charset=utf-8">\n')
f.write('<body>\n<ul>\n')
f.write('<style>table, th, td { border: 1px solid black; border-collapse: collapse; padding: 8px;}</style>\n')
f.write('<h1>Industrial Wireless Downloads</h1>')
introtext = 'This page contains links to all of the files currently available for download.  \
			 To download, click the link or right click and select Save-As.  Our downloadable \
			 files are stored on an external server which may cause some download difficulties. \
			 Utilities such as DownloadThemAll! are available to make downloading entire sections easier. \
			 <p>We have tried to make accessing the data as easy as possible given our resources. Please \
			 report broken links, missing data, or other requests to our project leader \
			 listed on our project home page.</p>'
f.write(introtext)

f.write("<h2><u>Table of Contents</u></h2>\n")

f.write('<h3>Summary files</h3>')
f.write('<table style="width:60%">')

f.write('<tr>')
href_figures = 'https://s3.amazonaws.com/nist-ics-wireless/measurements/active_meas/figures.zip'
f.write('<td>Channel estimation figures</td><td> <a href=%s>figures.zip</a></td>' % href_figures)
f.write('</tr>')

f.write('<tr>')
href_tr_cir = 'https://s3.amazonaws.com/nist-ics-wireless/measurements/active_meas/reduced_tap_delay_profiles.zip'
f.write('<td>13 tap average channel impulse response files</td><td>  <a href=%s>reduced_tap_delay_profiles.zip</a></td>' % href_tr_cir)
f.write('</tr>')

f.write('<tr>')
href_measst_cir = 'https://s3.amazonaws.com/nist-ics-wireless/measurements/active_meas/measurement_stats.zip'
f.write('<td>Channel estimates, metrics, and statistics</td><td>  <a href=%s>reduced_tap_delay_profiles.zip</a></td>' % href_measst_cir)
f.write('</tr>')

f.write('</table>')

f.write('<h3>Raw processed data files</h3>\n')
for ii in range(len(anchors)):
	f.write('<li><a href=#%s>%s</a><br></li>\n' % (anchors[ii], myGroupHeaders[ii]))
f.write('<br><hr>\n')
'''
f.write('<a name="%s"></a>' % 'reduced_tap_delay_profiles')
f.write('<u><h2>Reduced Tap Channel Impulse Response Files</h2></u>')
f.write('<p> Tap reduction algorithm is described in: </p><blockquote>C. Mehlfuhrer and M. Rupp, "Approximation and resampling of tapped delay \
		 line channel models with guaranteed channel properties," in 2008 IEEE \
		 International Conference on Acoustics, Speech and Signal Processing, 2008, \
		 no. 2, pp. 2869-2872.</blockquote>')'''
f.write('<u><h2>Raw Processed Data</h2></u>')
for ii in range(ndirs):
	print(ii)
	myDirPath = myDirPaths[ii]
	print("Working on %s ==> %s" % (myDirPath, myOutFile))
	myFiles = os.listdir(myDirPath)
	#print(myFiles)
	f.write('<a name="%s"></a>' % anchors[ii])
	f.write('<h3> %s </h3>' % myGroupHeaders[ii])
	f.writelines(['<a href="%s%s/%s">%s</a>...&nbsp\n' % (urlbase, myDirPath, f, f) for f in myFiles if f.endswith('.mat') or f.endswith('.csv') or f.endswith('.txt')])

f.write('<p>End of Downloads</p>\n')
f.write('</ul>\n</body>\n</htm>\n')
f.close()
