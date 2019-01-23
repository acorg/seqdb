# -*- Makefile -*-
# ----------------------------------------------------------------------

TARGETS = \
  $(SEQDB_LIB) \
  $(SEQDB_PY_LIB) \
  $(DIST)/seqdb-export-for-tree-maker \
  $(DIST)/seqdb-info \
  $(DIST)/seqdb-report-clade \
  $(DIST)/seqdb-report-dates \
  $(DIST)/seqdb-report-not-found-in-hidb \
  $(DIST)/seqdb-export-sequences-of-chart \
  $(DIST)/seqdb-chart-clades \
  $(DIST)/seqdb-clades-for-sera-in-chart \
  $(DIST)/seqdb-export-sequences-and-layout-of-chart \
  $(DIST)/seqdb-sequenced-antigens-in-chart \
  $(DIST)/seqdb-sequences-of-chart-for-ace-view \
  $(DIST)/seqdb-chart-antigen-to-seqid \
  $(DIST)/seqdb-strains-of-chart-in-clade \
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

all: install

CONFIGURE_PYTHON = 1
include $(ACMACSD_ROOT)/share/Makefile.config

LDLIBS = \
	$(AD_LIB)/$(call shared_lib_name,libacmacsbase,1,0) \
	$(AD_LIB)/$(call shared_lib_name,liblocationdb,1,0) \
	$(AD_LIB)/$(call shared_lib_name,libacmacschart,2,0) \
	$(AD_LIB)/$(call shared_lib_name,libhidb,5,0) \
	$(XZ_LIBS) $(PYTHON_LIBS) $(CXX_LIBS)

# ----------------------------------------------------------------------

install: install-headers $(TARGETS)
	$(call install_lib,$(SEQDB_LIB))
	$(call install_py_lib,$(SEQDB_PY_LIB))
	$(call symbolic_link_wildcard,$(DIST)/seqdb-*,$(AD_BIN))
	$(call symbolic_link_wildcard,$(abspath py)/*,$(AD_PY))
	$(call symbolic_link_wildcard,$(abspath bin)/seqdb-*,$(AD_BIN))
	$(call symbolic_link_wildcard,$(abspath bin)/fasta-*,$(AD_BIN))

test: install
	test/test
.PHONY: test

# ----------------------------------------------------------------------

$(SEQDB_LIB): $(patsubst %.cc,$(BUILD)/%.o,$(SEQDB_SOURCES)) | $(DIST) install-headers
	$(call echo_shared_lib,$@)
	$(call make_shared_lib,$(SEQDB_LIB_NAME),$(SEQDB_LIB_MAJOR),$(SEQDB_LIB_MINOR)) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(SEQDB_PY_LIB): $(patsubst %.cc,$(BUILD)/%.o,$(SEQDB_PY_SOURCES)) | $(DIST)
	$(call echo_shared_lib,$@)
	$(call make_shared_lib,$(SEQDB_PY_LIB_NAME),$(SEQDB_PY_LIB_MAJOR),$(SEQDB_PY_LIB_MINOR)) $(LDFLAGS) -o $@ $^ $(LDLIBS) $(PYTHON_LIBS)

$(DIST)/%: $(BUILD)/%.o | $(SEQDB_LIB)
	$(call echo_link_exe,$@)
	$(CXX) $(LDFLAGS) -o $@ $^ $(SEQDB_LIB) $(LDLIBS) $(AD_RPATH)

# ======================================================================
### Local Variables:
### eval: (if (fboundp 'eu-rename-buffer) (eu-rename-buffer))
### End:
