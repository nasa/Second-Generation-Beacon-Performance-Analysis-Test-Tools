NASA OPEN SOURCE AGREEMENT VERSION 1.3

THIS OPEN SOURCE AGREEMENT ("AGREEMENT") DEFINES THE RIGHTS OF USE,
REPRODUCTION, DISTRIBUTION, MODIFICATION AND REDISTRIBUTION OF CERTAIN
COMPUTER SOFTWARE ORIGINALLY RELEASED BY THE UNITED STATES GOVERNMENT
AS REPRESENTED BY THE GOVERNMENT AGENCY LISTED BELOW ("GOVERNMENT
AGENCY").  THE UNITED STATES GOVERNMENT, AS REPRESENTED BY GOVERNMENT
AGENCY, IS AN INTENDED THIRD-PARTY BENEFICIARY OF ALL SUBSEQUENT
DISTRIBUTIONS OR REDISTRIBUTIONS OF THE SUBJECT SOFTWARE.  ANYONE WHO
USES, REPRODUCES, DISTRIBUTES, MODIFIES OR REDISTRIBUTES THE SUBJECT
SOFTWARE, AS DEFINED HEREIN, OR ANY PART THEREOF, IS, BY THAT ACTION,
ACCEPTING IN FULL THE RESPONSIBILITIES AND OBLIGATIONS CONTAINED IN
THIS AGREEMENT.

Government Agency: NASA 
Government Agency Original Software Designation: GSC-18375-1
Government Agency Original Software Title: Second Generation Beacon Performance Analysis Test Tools
User Registration Requested.  Please Visit https://software.nasa.gov/
Government Agency Point of Contact for Original Software: David Watson, david.w.watson-1@nasa.gov




This is the release of the NASA T21 tools.  Create a development directory where you wish to save the files.
Extract the files into the chosen directory.  After extracting the files, run the top level file "doit.m", 
It will process an example file "testburst.mat" and produce the results: figures and excel measurements.

Once you have collected your own burst files of the format described in the document "T21 Signal Processing Scripts
Description_Release02032020.docx" you are ready to configure the scripts to process.

Open the file doit.m

modify the variable "xdir" and point it to the path where your collected files reside.

For example:

xdir = 'c:\mypath';

where files with the extension '*.mat' may be found.


If there is known attenuation in the receive path to the test equipment, you can account for this by modifying the 
variable "pcal".

For example, if a 20dB attenuator is placed between the beacon and test equipment, set pcal:

pcal=20;

Now when the power is measured it will reflect the EIRP of the beacon.

Once the script is modified with your path and power calibration inputs, save the script and run.

The scripts will read all of the mat files and process them.  Subdirectories will be created in the xdir path 
each of which contains the images of the measurements made on the burst.  In addition excel files will be generated
that contain the measurements and burst message payload.

Enjoy
Reese Bovard
Concentric Real Time, LLC

