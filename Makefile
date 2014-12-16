DEBUG?= -O2
PREFIX?=/usr/local/bin
CASMAP_ROOT?=/usr/local/share/casmap
CFLAGS?= $(DEBUG) -Wall -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -pthread  
LDFLAGS?= $(DEBUG) -Wall -lm -pthread  
BIN_PATH=$(CASMAP_ROOT)/bin
CC=gcc
BIN=cm_aln cm_ref cm_aln2pos cm_pos2sam cm_load cm_seq cm_sa cm_seq2sam cm_pos2sampe cm_sadbd
cm_ref-objs=seq.o refdb.o cm_ref.o
cm_sa-objs=cm_sa.o sadb.o alndb.o utils.o
cm_seq-objs=cm_seq.o seqdb.o
cm_load-objs=cm_load.o load.o utils.o 
cm_aln-objs=cm_aln.o aln.o seq.o seqdb.o alndb.o utils.o
cm_aln2pos-objs=seq.o sadb.o utils.o cm_aln2pos.o
cm_pos2sam-objs=seq.o refdb.o seqdb.o posdb.o pos2sam.o cm_pos2sam.o utils.o ksw.o global_aln.o  pos2sampe.o
cm_seq2sam-objs=seq.o refdb.o seqdb.o sadb.o aln.o aln2pos.o pos2sam.o pos2sampe.o cm_seq2sam.o global_aln.o  utils.o pipe.o  ksw.o
cm_pos2sampe-objs=seq.o refdb.o seqdb.o posdb.o cm_pos2sampe.o global_aln.o  pos2sampe.o utils.o ksw.o
cm_sadbd-objs=cm_sadbd.o sadb.o utils.o
CMD=casmap casmap-gui

%.pyc: %.py
	$(PYC) $^

all: $(BIN) 
#py-build $(CMD)

py-build:
	make -C py

target-build:
	make -C target

cm_aln: $(cm_aln-objs)
	$(CC) -o $@ $^ $(LDFLAGS) 

cm_ref: $(cm_ref-objs)
	$(CC) -o $@ $^ $(LDFLAGS) 

cm_aln2pos: $(cm_aln2pos-objs)
	$(CC) -o  $@ $^ $(LDFLAGS) 

cm_pos2sam: $(cm_pos2sam-objs)
	$(CC) -o $@ $^ $(LDFLAGS) 

cm_load: $(cm_load-objs)
	$(CC) -o $@ $^ $(LDFLAGS) 

cm_seq: $(cm_seq-objs)
	$(CC) -o $@ $^ $(LDFLAGS) 

cm_sa: $(cm_sa-objs)
	$(CC) -o $@ $^ $(LDFLAGS) 

cm_seq2sam: $(cm_seq2sam-objs)
	$(CC) -o $@ $^ $(LDFLAGS)

cm_pos2sampe: $(cm_pos2sampe-objs)
	$(CC) -o $@ $^  $(LDFLAGS) 

cm_sadbd: $(cm_sadbd-objs)
	$(CC) -o $@ $^ $(LDFLAGS)

casmap: casmap.sh
	echo "#!/bin/sh" > $@
	echo "CASMAP_ROOT=$(CASMAP_ROOT)" >> $@
	cat casmap.sh >> $@
	chmod 755 $@

casmap-gui: casmap-gui.sh
	echo "#!/bin/sh" > $@
	echo "CASMAP_ROOT=$(CASMAP_ROOT)" >> $@
	cat casmap-gui.sh >> $@
	chmod 755 $@

install: $(BIN)
	if [ ! -d $(CASMAP_ROOT) ];then mkdir $(CASMAP_ROOT);fi
	if [ ! -d $(BIN_PATH) ];then mkdir $(BIN_PATH);fi
	for a in $(BIN);do install $$a $(BIN_PATH)/$$a;done
	for a in $(CMD);do install $$a $(PREFIX)/$$a;done
	make -C py install

target-install: target-build
	make -C target install

uninstall:
	for a in $(CMD);do rm -f $(PREFIX)/$$a;done
	rm -rf $(CASMAP_ROOT)

clean:
##	make -C py clean
	rm -f $(CMD)
	rm -f $(BIN)
	rm -f *o

target-clean:
	make -C target clean

