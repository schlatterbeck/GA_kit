# To use this Makefile, get a copy of my SF Release Tools
# git clone git://git.code.sf.net/p/sfreleasetools/code sfreleasetools
# And point the environment variable RELEASETOOLS to the checkout

ifeq (,${RELEASETOOLS})
    RELEASETOOLS=../releasetools
endif
PKG=advanced-GA
PY=sga.py ecga.py hboa.py deceptive.py
README=README.rst
SRC=Makefile MANIFEST.in setup.py $(README) README.html $(PY)

VERSIONPY=Version.py
VERSION=$(VERSIONPY)
LASTRELEASE:=$(shell $(RELEASETOOLS)/lastrelease -n)

USERNAME=schlatterbeck
PROJECT=advanced-GA
PACKAGE=${PKG}

all: $(VERSION)

$(VERSION): $(SRC)

dist: all
	python setup.py sdist --formats=gztar,zip

clean:
	rm -rf default.css Version.py Version.pyc ${CLEAN}

include $(RELEASETOOLS)/Makefile-sf
