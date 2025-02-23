********************************************************
	notes on SITRINEO DAQ package
********************************************************

Contents:
---------
I.   RUNNING
II.  CONFIGURE
III. INSTALLING
IV.  DESCRIPTION
V.   VERSIONS
VI.  STRASBOURG SETUP



********************************************************
********************************************************
I.  RUNNING
********************************************************

You can only run the DAQ system if the previous INSTALL and CONFIGURE steps have been completed
and the DAQboard is reachable via ethernet from your computer.
That means you can ssh to the DAQboard and the initialization script knows your sensor configuration.

1/ Powering

Plug all cables:
-the cable with the banana plugs should be connected to the power supply (here TENMA 72-2540), the green plug in the + entry, the yellow in the - entry and the blue not connected
-the black cable should connect the DAQboard and the power outlet
-Ethernet should connect PC and DAQboard
Turn on the power supply, press the Off/On button (current ~ 1.5 A (for 4 sensors))
Wait a moment (10 seconds)
Press the red button on the DAQboard, and wait until the red lights blink 1 time 
Wait a moment (10 seconds)
Now everything is On

2/ Connect

Make sure your computer IP address has been changed to the DAQboard domain (usually 192.168.1.xxx).
See III. for explanation on the network configuration.
Ssh to the DAQboard as root (the IP address is written on the DAQboard): ssh root@192.168.1.12
If you get the error message : "No route to host", wait a bit longer until you retry
You can work from the root directory.

3/ Init and configure
Issue the following three instructions (see II. for explanations)
  source RESET_pulse_script
  daqSoC_udp -RunNumber 0 -NumEventsToRead 00 -triggerSW no -StepTriggerMonitor 50 -JtagInit work -IPnumber 0.0.0.0 -Delay 1
  source START_pulse_script
The consumption on the power supply should rise to ~ 2.2 A (for 4 sensors).
You should get among the messages something like (if 4 channels selected):
*************************************************************
******   Num Wrong Event Channel 0 =     0   Enabled   ******
******   Num Wrong Event Channel 1 =     0   Enabled   ******
******   Num Wrong Event Channel 2 =     0   Enabled   ******
******   Num Wrong Event Channel 3 =     0   Enabled   ******
*************************************************************


4/ Take data
Issue the following command with 3 arguments [RRRR], [EEEE], [MMMM]
  daqSoC_udp -RunNumber [RRRR] -NumEventsToRead [EEEE] -triggerSW yes -StepTriggerMonitor [MMMM] -JtagInit work -DataSaveLocal yes -Delay 1
where:
 RRRR = run number
 EEEE = number of events to take
 MMMM = drive the frequency of printouts, one line for monitoring written on stdout every MMMM events

It is always good to start with a small number of events for checking the read-out is fine, typically with EEEE~1000 and MMMM~50.
CHECK with no source that the number of words are low (typically 20 or lower) and exactly the same on both outputs of each sensor
Note that the first line (event number 0) contains always crap.
NumEvent-counter-Busy reg *  ch0   - ch0   - ch1   - ch1   -  ch2   - ch2   - ch3  - ch3    * ch0 -ch0 -ch1 -ch1 -ch2 -ch2 -ch3 -ch3
NumEvent-counter-Busy reg *  DO0   - DO1   - DO0   - DO1   -  DO0   - DO1   - DO0  - DO1    * DO0 -DO1 -DO0 -DO1 -DO0 -DO1 -DO0 -DO1
       0       1   ff0000 *     504    1557    1085      44     971     971     494    1069 *    1    1    1    1 5068 5068    1    1
     150     151   ff0000 *     504    1557    1085     899    1637    1637    1103    1678 *    1    1  677    1   21   21   19   19
If you detect troubles like in the example above, look at II. CONFIGURING Phase shift

Now you can take more data...enjoy!

********************************************************
********************************************************
II. CONFIGURING
********************************************************

** Sensor choice
daqSoC_init_script
At the bottom of the file:
accessHPS setbits   ff200000  40000    10000 1 > /dev/null  # enable ch 0
accessHPS setbits   ff200000  40000    20000 1 > /dev/null  # enable ch 1
accessHPS setbits   ff200000  40000    40000 1 > /dev/null  # enable ch 2
accessHPS setbits   ff200000  40000    80000 1 > /dev/null  # enable ch 3


** JTAG
daqSoC_Jtag_init_script_work
Example for channel 0 with 4 threshold values (aka Vref1), 198, 157, 183 and 183:
#            Channel 0
jtagProgrammingSoC_v3 -action WriteSubRegister -DRbitStreamFile JtagBitsStreamFile_work.dat -Register BIAS_GEN -SubRegister v_disc_ref1A -Format Dec -dataToWriteDec 198 > /dev/null
jtagProgrammingSoC_v3 -action WriteSubRegister -DRbitStreamFile JtagBitsStreamFile_work.dat -Register BIAS_GEN -SubRegister v_disc_ref1B -Format Dec -dataToWriteDec 157 > /dev/null
jtagProgrammingSoC_v3 -action WriteSubRegister -DRbitStreamFile JtagBitsStreamFile_work.dat -Register BIAS_GEN -SubRegister v_disc_ref1C -Format Dec -dataToWriteDec 183 > /dev/null
jtagProgrammingSoC_v3 -action WriteSubRegister -DRbitStreamFile JtagBitsStreamFile_work.dat -Register BIAS_GEN -SubRegister v_disc_ref1D -Format Dec -dataToWriteDec 183 > /dev/null

Rule of thumb to set the threshold: 4 digital units = 1 sigma (of noise)
Usually, the initial setting is for threshold = 8*sigma + pedestal = 32 digits + middle_point
Maximal value is 255


** Phase shift
 - sensor 2 (channels 3-4) phase-shift 2
	PHASESHIFT_pulse_script 3 2
	PHASESHIFT_pulse_script 4 2


********************************************************
********************************************************
III. INSTALLING
********************************************************

/ Network configuration
Configure ethernet address on
sudo ifconfig eth0 192.168.1.100 netmask 255.255.255.0


********************************************************
********************************************************
IV. DESCRIPTION
********************************************************




********************************************************
********************************************************
V. VERSION
********************************************************
-IPnumberUDPmon 0.0.0.0 for old data format
-DataSaveLocal yes USE -IPnumberUDP for new data format


********************************************************
********************************************************
VI.  STRASBOURG SETUP
********************************************************

* TENMA power supply, memory M1 with Sitrineo configuration (voltage 5 V, current limit 2.5 A)

* The network board on the local machine is eno2 and is registered as 192.168.1.2
sudo ifconfig eno2 192.168.1.2 netmask 255.255.255.0

* To open a connection with the DAQboard
ssh root@192.168.1.12

There is an alias "sitrineo" in the Terminal App to connect directly to the DAQboard

* Shortcuts on the DAQboard
sitrireset => RESET
sitriconf => JTAG configuration
sitristart => start the clockconf
sitridaq => takes a short run (not saving data)

* Taking data stored on computer 
Go to daq directory and use command (not working for the moment)
./m2DAQ.sh [runNumber] [nbOfEvents]

other possibility, take data on the Sitrineo board and then on PC under ~/data directory
scp -r root@192.168.1.12:DataStore/run0001 1
(for run number 1)

* Tips to send files from daqBoard to PC
scp root@192.168.1.12:myfile .

* Tips to send files from PC to daqBoard
scp myfile root@192.168.1.12:.

** Using TAF
* Once per session
source Scripts/thistaf.sh 

* Launch analysis
taf -run 1 -cfg ./config/sitrineo-m2.cfg -guiw
and click on the menu... Cumulate Raw Data 2D

* TAF handy Raw & Hit analysis commands
gTAF->GetRaw()->DisplayRawData2D()
gTAF->GetRaw()->DisplayCumulatedHits2D(nEvents) // ex. (10000)

* TAF handy Alignment commands:
gTAF->SetAlignStatus(0) // Set 1st plane as reference plane
gTAF->SetAlignStatus(1) // Set 2nd plane as reference plane
gTAF->AlignTracker(Area in um, 0, nEvents) // ex. (1000,0,20000)



