# About
The LimeSDR does not interface with GNU Radio as nicely as we had hoped. I imagine that at some point in the future things will smooth out but until then this is what might help.

Following are instruction sets for:
- Linux
- macOS
- Windows (not yet)

Take everything with a grain of salt. I'm a grad student, I'm not swimming in money and spare time so I rigged up something that works for me and I'm willing to share that solution here.

## Background reading
*stuff you should read about before launching into the install so I can just say things like "Soapy is good" and you don't think I'm a moron.*

# Linux
This was tested with Ubuntu 16.04 and LimeSDR-USB with:
- Firmware: 3
- Hardware: 4
- Protocol: 1
- Gateware: 2
- Gateware rev.: 8

Using the apt-get versions of gnuradio/soapysdr/limesuite/gr-osmosdr caused failure problems for me, and there is conflicting info across the Discourse forum about how to get things working. I ended up using [this post](https://discourse.myriadrf.org/t/latest-aug-2017-limesdr-and-gnuradio/1563) as my primary guide.

I started with nothing installed at all. I wiped any sdr related PPAs from my system, and everything on my system that included the strings "soapy", "lime", "uhd", "gnuradio", "grc", "bladerf", "rtl-sdr", "gr-osmosdr", "myriadrf", "osmo", "hackrf", "umtrx", "miri", "rfspace", "ettus" or "pothos".

## Package manager
Pybombs is what I used to get everything installed and interfacing happily.
```
sudo apt install pip python-apt
sudo pip install --upgrade pip
sudo pip install pybombs
```

PyBombs needs recipe files:
```
pybombs recipes add gr-recipes git+https://github.com/gnuradio/gr-recipes.git
```

I made a directory to house GNU Radio in my users home folder and started the install (it took a while).
```
mkdir prefix
pybombs prefix init -a default prefix/default/ -R gnuradio-default
```

If anything fails along the way, you can restart the process with just:
```
pybombs install gnuradio
```

Once that has been done successfully, get the fun stuff for the Limes to use.
```
pybombs install soapysdr limesuite gr-osmosdr
```

Use the handily included environment configuration script to set a few environment variables:
```
source /prefix/default/setup-env.sh
```

I had to change some permissions to allow the LimeUtil and SoapySDRUtil programs to read/write to USB with:
```
cd prefix/default/src/limesuite/udev-rules/
sudo ./install.sh
```

The I tested the LimeSDR by running each of the following commands and making sure that everything worked fine.
```
LimeUtil --find
SoapySDRUtil --find
```


# macOS
This was tested with macOS High Sierra (10.13.1) and LimeSDR-USB with:
- Firmware: 3
- Hardware: 4
- Protocol: 1
- Gateware: 2
- Gateware rev.: 10


## Xcode Command Line Tools
Pretty much everything requires this. You may already have it installed. If you run:
```
xcode-select --install
```
It'll either launch an install dialog or say "go away this is already installed". If it launches the dialog choose "Install" not "Get Xcode". This install will take a while (think hours depending on your connection).

# _Archived:_
In the past a large installation process was necessary but if you skip down to the [GNURadio installation instructions](#gnu-radio) there is now a way to avoid the majority of it.
> ## Package manager
> I've classically used MacPorts for everything but for this I went with Homebrew. Installing GNU Radio with Macports is supported and easy, but installing everything else is only well supported via brew. So take your pick. I opted for annoying gnuradio install and easy everything else install.

> **Do not have both MacPorts and HomeBrew installed at once** - It'll just suck, probably. To uninstall MacPorts go to [their guide](https://guide.macports.org/chunked/installing.macports.uninstalling.html), they update it regularly so I don't want to copy out the commands here.


> ### Install homebrew
> Visit [their site](https://brew.sh/) where there will be up-to-date notes if you're on a weird or newer OS. The probably permanently useful install can be executed with:
> ```
> /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
> ```

> ## LimeSuite and Soapy
> This is the stuff that lets the board communicate with your OS and vice versa. They support brew natively which is rad and their install guide can be followed [here](https://wiki.myriadrf.org/Lime_Suite).

> For convenience:
> ```
> brew tap pothosware/homebrew-pothos
> brew update
> brew install limesuite
> ```

> To test that this worked, plug your LimeSDR into the computer and run each of:
> ```
> LimeUtil --find
> ```
> ```
> SoapySDRUtil --find
> ```
> They should each return info about your board.

## GNU Radio
[This handy repo](https://github.com/cfriedt/gnuradio-for-mac-without-macports) allows you to install GNURadio as a `.app`, and as of July it includes support for the LimeSDR!

It does require python 2.7 and although you probably definitely already have it installed on your system I went ahead and installed it again to make sure everything was up to date.

### Testing
At this point in theory everything is in place. Launch GNU Radio, plug in your Lime, and try out one of our test scripts located in <https://github.com/nsbruce/UViip/tree/master/GNUradio/LimeSDR>.
