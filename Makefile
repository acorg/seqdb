# -*- Makefile -*-
# Eugene Skepner 2017
# ----------------------------------------------------------------------

MAKEFLAGS = -w

# ----------------------------------------------------------------------

TARGETS = \
	$(SEQDB_LIB) \
	$(SEQDB_PY_LIB) \
	$(DIST)/seqdb-info \
	$(DIST)/seqdb-report-clade \
	$(DIST)/seqdb-report-dates \
	$(DIST)/seqdb-report-not-found-in-hidb \
	$(DIST)/seqdb-export-sequences-of-chart \
	$(DIST)/seqdb-export-sequences-and-layout-of-chart \
	$(DIST)/seqdb-update-clades \
	$(DIST)/seqdb-list-strains-in-the-clade \
	$(DIST)/seqdb-list-strains-having-aa-at \
	$(DIST)/seqdb-compare-sequences

SEQDB_SOURCES = seqdb.cc seqdb-export.cc seqdb-import.cc seqdb-hidb.cc amino-acids.cc clades.cc insertions_deletions.cc
SEQDB_PY_SOURCES = $(SEQDB_SOURCES) py.cc

SEQDB_LIB_MAJOR = 2
SEQDB_LIB_MINOR = 0
SEQDB_LIB_NAME = libseqdb
SEQDB_LIB = $(DIST)/$(call shared_lib_name,$(SEQDB_LIB_NAME),$(SEQDB_LIB_MAJOR),$(SEQDB_LIB_MINOR))

SEQDB_PY_LIB_MAJOR = 1
SEQDB_PY_LIB_MINOR = 0
SEQDB_PY_LIB_NAME = seqdb_backend
SEQDB_PY_LIB = $(DIST)/$(SEQDB_PY_LIB_NAME)$(PYTHON_MODULE_SUFFIX)

# ----------------------------------------------------------------------

include $(ACMACSD_ROOT)/share/makefiles/Makefile.g++
include $(ACMACSD_ROOT)/share/makefiles/Makefile.python
include $(ACMACSD_ROOT)/share/makefiles/Makefile.dist-build.vars

CXXFLAGS = -MMD -g $(OPTIMIZATION) $(PROFILE) -fPIC -std=$(STD) $(WARNINGS) -I$(AD_INCLUDE) $(PKG_INCLUDES)
LDFLAGS = $(OPTIMIZATION) $(PROFILE)
LDLIBS = \
	$(AD_LIB)/$(call shared_lib_name,libacmacsbase,1,0) \
	$(AD_LIB)/$(call shared_lib_name,liblocationdb,1,0) \
	$(AD_LIB)/$(call shared_lib_name,libacmacschart,2,0) \
	$(AD_LIB)/$(call shared_lib_name,libhidb,5,0) \
	$(shell pkg-config --libs liblzma) \
	$(shell $(PYTHON_CONFIG) --ldflags | sed -E 's/-Wl,-stack_size,[0-9]+//') \
	$(CXX_LIB)

PKG_INCLUDES = $(shell pkg-config --cflags liblzma) $(PYTHON_INCLUDES)

# ----------------------------------------------------------------------

all: check-acmacsd-root install-headers $(TARGETS)

install: check-acmacsd-root install-headers $(TARGETS)
	$(call install_lib,$(SEQDB_LIB))
	$(call install_py_lib,$(SEQDB_PY_LIB))
	ln -sf $(abspath dist)/seqdb-* $(AD_BIN)
	ln -sf $(abspath py)/* $(AD_PY)
	ln -sf $(abspath bin)/seqdb-* $(AD_BIN)

test: install
	test/test

# ----------------------------------------------------------------------

-include $(BUILD)/*.d
include $(ACMACSD_ROOT)/share/makefiles/Makefile.dist-build.rules
include $(ACMACSD_ROOT)/share/makefiles/Makefile.rtags

# ----------------------------------------------------------------------

$(SEQDB_LIB): $(patsubst %.cc,$(BUILD)/%.o,$(SEQDB_SOURCES)) | $(DIST) install-headers
	@printf "%-16s %s\n" "SHARED" $@
	@$(call make_shared,$(SEQDB_LIB_NAME),$(SEQDB_LIB_MAJOR),$(SEQDB_LIB_MINOR)) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(SEQDB_PY_LIB): $(patsubst %.cc,$(BUILD)/%.o,$(SEQDB_PY_SOURCES)) | $(DIST)
	@printf "%-16s %s\n" "SHARED" $@
	@$(call make_shared,$(SEQDB_PY_LIB_NAME),$(SEQDB_PY_LIB_MAJOR),$(SEQDB_PY_LIB_MINOR)) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(PYTHON_LDLIBS)

$(DIST)/%: $(BUILD)/%.o | $(SEQDB_LIB)
	@printf "%-16s %s\n" "LINK" $@
	@$(CXX) $(LDFLAGS) -o $@ $^ $(SEQDB_LIB) $(LDLIBS) $(AD_RPATH)

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
