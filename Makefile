# -*- Makefile -*-
# Eugene Skepner 2017
# ----------------------------------------------------------------------

MAKEFLAGS = -w

# ----------------------------------------------------------------------

SEQDB_SOURCES = seqdb.cc seqdb-export.cc seqdb-import.cc seqdb-hidb.cc amino-acids.cc clades.cc insertions_deletions.cc
SEQDB_PY_SOURCES = $(SEQDB_SOURCES) py.cc
SEQDB_REPORT_CLADE_SRC = seqdb-report-clade.cc
SEQDB_REPORT_DATES_SRC = seqdb-report-dates.cc

# ----------------------------------------------------------------------

include $(ACMACSD_ROOT)/share/makefiles/Makefile.g++
include $(ACMACSD_ROOT)/share/makefiles/Makefile.dist-build.vars

PYTHON_VERSION = $(shell python3 -c 'import sys; print("{0.major}.{0.minor}".format(sys.version_info))')
PYTHON_CONFIG = python$(PYTHON_VERSION)-config
PYTHON_MODULE_SUFFIX = $(shell $(PYTHON_CONFIG) --extension-suffix)

SEQDB_LIB = $(DIST)/libseqdb.so

CXXFLAGS = -MMD -g $(OPTIMIZATION) $(PROFILE) -fPIC -std=$(STD) $(WARNINGS) -I$(AD_INCLUDE) $(PKG_INCLUDES)
LDFLAGS = $(OPTIMIZATION) $(PROFILE)
SEQDB_LDLIBS = -L$(AD_LIB) -lacmacsbase -lacmacschart -llocationdb -lhidb -lacmacsbase $(shell pkg-config --libs liblzma) $(shell $(PYTHON_CONFIG) --ldflags | sed -E 's/-Wl,-stack_size,[0-9]+//')

PKG_INCLUDES = $(shell pkg-config --cflags liblzma) $(shell $(PYTHON_CONFIG) --includes)

# ----------------------------------------------------------------------

all: check-acmacsd-root $(DIST)/seqdb_backend$(PYTHON_MODULE_SUFFIX) $(SEQDB_LIB) $(DIST)/seqdb-report-clade $(DIST)/seqdb-report-dates

install: check-acmacsd-root install-headers $(DIST)/seqdb_backend$(PYTHON_MODULE_SUFFIX) $(SEQDB_LIB) $(DIST)/seqdb-report-clade $(DIST)/seqdb-report-dates
	ln -sf $(DIST)/seqdb_backend$(PYTHON_MODULE_SUFFIX) $(AD_PY)
	ln -sf $(DIST)/seqdb-report-* $(AD_BIN)
	ln -sf $(abspath py)/* $(AD_PY)
	ln -sf $(abspath bin)/seqdb-* $(AD_BIN)

install-libseqdb: $(SEQDB_LIB)
	$(call install_lib,$^)

test: install
	test/test

# ----------------------------------------------------------------------

-include $(BUILD)/*.d
include $(ACMACSD_ROOT)/share/makefiles/Makefile.dist-build.rules
include $(ACMACSD_ROOT)/share/makefiles/Makefile.rtags

# ----------------------------------------------------------------------

$(DIST)/seqdb_backend$(PYTHON_MODULE_SUFFIX): $(patsubst %.cc,$(BUILD)/%.o,$(SEQDB_PY_SOURCES)) | $(DIST) install-headers
	@echo "SHARED     " $@ # '<--' $^
	@$(CXX) -shared $(LDFLAGS) -o $@ $^ $(SEQDB_LDLIBS)
	@#strip $@

$(SEQDB_LIB): $(patsubst %.cc,$(BUILD)/%.o,$(SEQDB_SOURCES)) | $(DIST) $(LOCATION_DB_LIB) install-headers
	@echo "SHARED     " $@ # '<--' $^
	@$(CXX) -shared $(LDFLAGS) -o $@ $^ $(SEQDB_LDLIBS)

$(DIST)/seqdb-report-clade: $(patsubst %.cc,$(BUILD)/%.o,$(SEQDB_REPORT_CLADE_SRC)) | $(DIST) install-libseqdb install-headers
	@echo "LINK       " $@
	@$(CXX) $(LDFLAGS) -o $@ $^ -lseqdb $(SEQDB_LDLIBS)

$(DIST)/seqdb-report-dates: $(patsubst %.cc,$(BUILD)/%.o,$(SEQDB_REPORT_DATES_SRC)) | $(DIST) install-libseqdb install-headers
	@echo "LINK       " $@
	@$(CXX) $(LDFLAGS) -o $@ $^ -lseqdb $(SEQDB_LDLIBS)

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
