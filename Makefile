# -*- Makefile -*-
# Eugene Skepner 2016
# ----------------------------------------------------------------------

MAKEFLAGS = -w

# ----------------------------------------------------------------------

SEQDB_SOURCES = seqdb.cc seqdb-export.cc seqdb-import.cc seqdb-hidb.cc amino-acids.cc clades.cc insertions_deletions.cc
SEQDB_PY_SOURCES = $(SEQDB_SOURCES) py.cc
SEQDB_REPORT_CLADE_SRC = seqdb-report-clade.cc

# ----------------------------------------------------------------------

TARGET_ROOT=$(shell if [ -f /Volumes/rdisk/ramdisk-id ]; then echo /Volumes/rdisk/AD; else echo $(ACMACSD_ROOT); fi)
include $(TARGET_ROOT)/share/Makefile.g++
include $(TARGET_ROOT)/share/Makefile.dist-build.vars

PYTHON_VERSION = $(shell python3 -c 'import sys; print("{0.major}.{0.minor}".format(sys.version_info))')
PYTHON_CONFIG = python$(PYTHON_VERSION)-config
PYTHON_MODULE_SUFFIX = $(shell $(PYTHON_CONFIG) --extension-suffix)

SEQDB_LIB = $(DIST)/libseqdb.so

# -fvisibility=hidden and -flto make resulting lib smaller (pybind11) but linking is much slower
OPTIMIZATION = -O3 #-fvisibility=hidden -flto
PROFILE = # -pg
CXXFLAGS = -MMD -g $(OPTIMIZATION) $(PROFILE) -fPIC -std=$(STD) $(WEVERYTHING) $(WARNINGS) -I$(AD_INCLUDE) $(PKG_INCLUDES)
LDFLAGS = $(OPTIMIZATION) $(PROFILE)
SEQDB_LDLIBS = -L$(AD_LIB) -lacmacsbase -lacmacschart -llocationdb -lhidb -lacmacsbase -lboost_filesystem -lboost_system $(shell pkg-config --libs liblzma) $(shell $(PYTHON_CONFIG) --ldflags | sed -E 's/-Wl,-stack_size,[0-9]+//')

PKG_INCLUDES = $(shell pkg-config --cflags liblzma) $(shell $(PYTHON_CONFIG) --includes)

# ----------------------------------------------------------------------

all: check-acmacsd-root $(DIST)/seqdb_backend$(PYTHON_MODULE_SUFFIX) $(SEQDB_LIB) $(DIST)/seqdb-report-clade

install: check-acmacsd-root install-headers $(DIST)/seqdb_backend$(PYTHON_MODULE_SUFFIX) $(SEQDB_LIB) $(DIST)/seqdb-report-clade
	ln -sf $(SEQDB_LIB) $(AD_LIB)
	if [ $$(uname) = "Darwin" ]; then /usr/bin/install_name_tool -id $(AD_LIB)/$(notdir $(SEQDB_LIB)) $(AD_LIB)/$(notdir $(SEQDB_LIB)); fi
	ln -sf $(DIST)/seqdb_backend$(PYTHON_MODULE_SUFFIX) $(AD_PY)
	ln -sf $(DIST)/seqdb-report-clade $(AD_BIN)
	ln -sf $(abspath py)/* $(AD_PY)
	ln -sf $(abspath bin)/seqdb-* $(AD_BIN)

install-headers:
	if [ ! -d $(AD_INCLUDE)/seqdb ]; then mkdir $(AD_INCLUDE)/seqdb; fi
	ln -sf $(abspath cc)/*.hh $(AD_INCLUDE)/seqdb

test: install
	test/test

# ----------------------------------------------------------------------

-include $(BUILD)/*.d

include $(AD_SHARE)/Makefile.rtags

# ----------------------------------------------------------------------

$(DIST)/seqdb_backend$(PYTHON_MODULE_SUFFIX): $(patsubst %.cc,$(BUILD)/%.o,$(SEQDB_PY_SOURCES)) | $(DIST)
	@echo "SHARED     " $@ # '<--' $^
	@$(CXX) -shared $(LDFLAGS) -o $@ $^ $(SEQDB_LDLIBS)
	@#strip $@

$(SEQDB_LIB): $(patsubst %.cc,$(BUILD)/%.o,$(SEQDB_SOURCES)) | $(DIST) $(LOCATION_DB_LIB)
	@echo "SHARED     " $@ # '<--' $^
	@$(CXX) -shared $(LDFLAGS) -o $@ $^ $(SEQDB_LDLIBS)

$(DIST)/seqdb-report-clade: $(patsubst %.cc,$(BUILD)/%.o,$(SEQDB_REPORT_CLADE_SRC)) | $(DIST) $(SEQDB_LIB)
	@echo "LINK       " $@
	@$(CXX) $(LDFLAGS) -o $@ $^ -lseqdb $(SEQDB_LDLIBS)

# ----------------------------------------------------------------------

$(BUILD)/%.o: cc/%.cc | $(BUILD) install-headers
	@echo $(CXX_NAME) $<
	@$(CXX) $(CXXFLAGS) -c -o $@ $<

# ----------------------------------------------------------------------

include $(AD_SHARE)/Makefile.dist-build.rules

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
