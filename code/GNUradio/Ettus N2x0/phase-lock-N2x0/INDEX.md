# Test 1
One USRP N210 used from *uru* to TX. Coax from tx port to BNC splitter. Two cables of equal length going to the two RX USRPs. They are run from *mistic-elm*. Each USRP connected via distinct ethernet cables to a switch.
 - tx1-GR.png - GR flowgraph for transmitting. Simple signal source to USRP Sink
 - tx1-time.png - Shows real and imag output from signal source in time domain
 - rx1-GR.png - GR flowgraph for receiving. Only __real__ part of signals are kept (removed imaginary to make it look cleaner and since we only care about phase right now). Their amplitude is normalized using a feed-forward loop in the AGC block.
 - rx1-time.png - Real signals from each antenna in time domain (normalized magnitude)

# Test 2
One USRP N210 used from *uru* to TX. Coax from tx port to BNC splitter. Two cables of equal length going to the two RX USRPs. They are run from *mistic-elm*. Each USRP connected via distinct ethernet cables to a switch.
- tx1-GR.png - No change
- tx2-time.png - Difference frequency used
- rx2-GR.png - Flow graph outputs phase difference.
- rx2-phase-difference.png - Obviously not phase locked.

# Test 3
One USRP N210 used from *uru* to TX. Coax from tx port to BNC splitter. Two cables of equal length going to the two RX USRPs. They are run from *mistic-elm*. One USRP connected via ethernet cable and the other MIMOd in.
 - tx1-GR.png - No change
 - tx2-time.png - No change
 - rx3-GR.png - Uses MIMO
 - rx3-phase-difference.png - Still not phase locked.

# Test 4
One USRP N210 used from *uru* to TX. Coax from tx port to BNC splitter. Two cables of equal length going to the two RX USRPs. They are run from *mistic-elm*. One USRP connected via ethernet cable and the other MIMOd in.
- tx1-GR.png - No change
- tx2-time.png - No change
- rx3-GR.png - Sync using PC Clock
- rx4-phase-difference.png - Still not phase locked, but somehow symmetrical.

# Test 5
Waveform generator used to Tx. Coax from tx port to BNC splitter. Two cables of equal length but unequal everything()?) going to the two RX USRPs. They are run from *mistic-elm*. Each USRP ethernetted into the switch. 10MHz out from Flex to CLK in on RX USRPs.
- tx at 8.1MHz
- rx5-GR.png - Everything references external
- rx5-phase-difference.png - Change in phase seems constant.

# Test 6
Waveform generator used to Tx. Coax from tx port to BNC splitter. Two equal cables going to the two RX USRPs. They are run from *mistic-elm*. USRPs MIMOd. 1 USRP ethernetted in.
- tx at 30kHz
- Although the signals look constant frequency there is some base lower frequency. The entire rx signal sink drifts up and down.
- rx6-GR.png - Everything references MIMO
- rx6-phase-difference.png - Phase constant with random spikes.

# Test 7
Waveform generator used to Tx. Coax from tx port to BNC splitter. Two equal cables going to the two RX USRPs. They are run from *mistic-elm*. Each USRP ethernetted in.
- tx at 30kHz
- Although the signals look constant frequency there is some base lower frequency. The entire rx signal sink drifts up and down.
- rx7-GR.png - Everything default
- rx7-phase-difference.png - Phase un-constant.

# Test 8
Waveform generator used to Tx. Coax from tx port to BNC splitter. Two _unequal_ cables going to the two RX USRPs. They are run from *mistic-elm*. Each USRP ethernetted in.
- tx at 30kHz
- Although the signals look constant frequency there is some base lower frequency. The entire rx signal sink drifts up and down.
- rx7-GR.png - Everything default
- rx8-phase-difference.png - Phase un-constant.

# Test 9
Waveform generator used to Tx. Coax from tx port to BNC splitter. Two _unequal_ cables going to the two RX USRPs. They are run from *mistic-elm*. Each USRP ethernetted in. Two equal cables going from Flex CLK out to BNC splitter, to each USRPs CLK IN.
- tx at 30kHz
- rx9-GR.png - Clock source set to "external"
- rx9-phase-difference.png - Phase un-constant.

# Test 10
Waveform generator used to Tx. Coax from tx port to BNC splitter. Two _unequal_ cables going to the two RX USRPs. They are run from *mistic-elm*. Each USRP ethernetted in. Two equal cables going from Flex CLK out to BNC splitter, to each USRPs CLK IN.
- tx at 30kHz
- rx10-GR.png - Clock source set to "external", sync set to "PC Clock"
- rx10-phase-difference.png - Phase un-constant.

# Test 11
Waveform generator used to Tx. Coax from tx port to 2 BNC splitters. 4 unequal cables going two to each USRP. They are run from *mistic-elm*. One USRP ethernetted in. USRPs mimoed
- tx at 100kHz
- rx11-GR.png - need to upload
- rx11-phase-difference.png - Phase constant all four ways!
