# -*- Makefile -*-
# Eugene Skepner 2016
# ----------------------------------------------------------------------

MAKEFLAGS = -w

# ----------------------------------------------------------------------

SEQDB_SOURCES = seqdb.cc seqdb-export.cc seqdb-import.cc seqdb-hidb.cc amino-acids.cc clades.cc
SEQDB_PY_SOURCES = $(SEQDB_SOURCES) seqdb-py.cc

# ----------------------------------------------------------------------

CLANG = $(shell if g++ --version 2>&1 | grep -i llvm >/dev/null; then echo Y; else echo N; fi)
ifeq ($(CLANG),Y)
  WEVERYTHING = -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded
  WARNINGS = -Wno-weak-vtables # -Wno-padded
  STD = c++14
else
  WEVERYTHING = -Wall -Wextra
  WARNINGS =
  STD = c++14
endif

PYTHON_VERSION = $(shell python3 -c 'import sys; print("{0.major}.{0.minor}".format(sys.version_info))')
PYTHON_CONFIG = python$(PYTHON_VERSION)-config
PYTHON_MODULE_SUFFIX = $(shell $(PYTHON_CONFIG) --extension-suffix)

LIB_DIR = $(ACMACSD_ROOT)/lib
SEQDB_LIB = $(DIST)/libseqdb.so

# -fvisibility=hidden and -flto make resulting lib smaller (pybind11) but linking is much slower
OPTIMIZATION = -O3 #-fvisibility=hidden -flto
PROFILE = # -pg
CXXFLAGS = -MMD -g $(OPTIMIZATION) $(PROFILE) -fPIC -std=$(STD) $(WEVERYTHING) $(WARNINGS) -I$(BUILD)/include -I$(ACMACSD_ROOT)/include $(PKG_INCLUDES)
LDFLAGS = $(OPTIMIZATION) $(PROFILE)
SEQDB_LDLIBS = -L$(LIB_DIR) -llocationdb -lhidb -lacmacsbase -lboost_filesystem -lboost_system $$(pkg-config --libs liblzma) $$($(PYTHON_CONFIG) --ldflags | sed -E 's/-Wl,-stack_size,[0-9]+//')

PKG_INCLUDES = $$(pkg-config --cflags liblzma) $$($(PYTHON_CONFIG) --includes)

# ----------------------------------------------------------------------

BUILD = build
DIST = $(abspath dist)

all: check-acmacsd-root $(DIST)/seqdb_backend$(PYTHON_MODULE_SUFFIX) $(SEQDB_LIB)

install: check-acmacsd-root install-headers $(DIST)/seqdb_backend$(PYTHON_MODULE_SUFFIX) $(SEQDB_LIB)
	ln -sf $(SEQDB_LIB) $(ACMACSD_ROOT)/lib
	if [ $$(uname) = "Darwin" ]; then /usr/bin/install_name_tool -id $(ACMACSD_ROOT)/lib/$(notdir $(SEQDB_LIB)) $(ACMACSD_ROOT)/lib/$(notdir $(SEQDB_LIB)); fi
	ln -sf $(DIST)/seqdb_backend$(PYTHON_MODULE_SUFFIX) $(ACMACSD_ROOT)/py
	ln -sf $(abspath py)/* $(ACMACSD_ROOT)/py
	ln -sf $(abspath bin)/seqdb-* $(ACMACSD_ROOT)/bin

install-headers:
	if [ ! -d $(ACMACSD_ROOT)/include/seqdb ]; then mkdir $(ACMACSD_ROOT)/include/seqdb; fi
	ln -sf $(abspath cc)/*.hh $(ACMACSD_ROOT)/include/seqdb

test: install
	test/test

# ----------------------------------------------------------------------

-include $(BUILD)/*.d

# ----------------------------------------------------------------------

$(DIST)/seqdb_backend$(PYTHON_MODULE_SUFFIX): $(patsubst %.cc,$(BUILD)/%.o,$(SEQDB_PY_SOURCES)) | $(DIST)
	g++ -shared $(LDFLAGS) -o $@ $^ $(SEQDB_LDLIBS)
	@#strip $@

$(SEQDB_LIB): $(patsubst %.cc,$(BUILD)/%.o,$(SEQDB_SOURCES)) | $(DIST) $(LOCATION_DB_LIB)
	g++ -shared $(LDFLAGS) -o $@ $^ $(SEQDB_LDLIBS)

clean:
	rm -rf $(DIST) $(BUILD)/*.o $(BUILD)/*.d

distclean: clean
	rm -rf $(BUILD)

# ----------------------------------------------------------------------

$(BUILD)/%.o: cc/%.cc | $(BUILD) install-headers
	@echo $<
	@g++ $(CXXFLAGS) -c -o $@ $<

# ----------------------------------------------------------------------

check-acmacsd-root:
ifndef ACMACSD_ROOT
	$(error ACMACSD_ROOT is not set)
endif

$(DIST):
	mkdir -p $(DIST)

$(BUILD):
	mkdir -p $(BUILD)

.PHONY: check-acmacsd-root

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
