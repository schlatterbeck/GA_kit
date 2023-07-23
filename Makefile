# To use this Makefile, get a copy of my SF Release Tools
# git clone git://git.code.sf.net/p/sfreleasetools/code sfreleasetools
# And point the environment variable RELEASETOOLS to the checkout

ifeq (,${RELEASETOOLS})
    RELEASETOOLS=../releasetools
endif
PKG=GA_kit
PY=$(PKG)/sga.py $(PKG)/ecga.py $(PKG)/hboa.py $(PKG)/deceptive.py
README=README.rst
SRC=Makefile MANIFEST.in setup.py $(README) README.html $(PY)

VERSIONPY=$(PKG)/Version.py
VERSION=$(VERSIONPY)
LASTRELEASE:=$(shell $(RELEASETOOLS)/lastrelease -n)

USERNAME=schlatterbeck
PROJECT=GA_kit
PACKAGE=${PKG}

all: $(VERSION)

$(VERSION): $(SRC)

clean:
	rm -rf default.css Version.py Version.pyc ${CLEAN}

include $(RELEASETOOLS)/Makefile-pyrelease
