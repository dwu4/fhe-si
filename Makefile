CC = g++

CFLAGS = -g -O2 -Wall -std=c++0x -I/usr/local/include -Wno-delete-non-virtual-dtor

BUILD = build
TESTS = tests

LDPATH = -L/usr/lib/ -L$(CURDIR)/$(BUILD)
LDLIBS = -lntl -lm

SRC =  PlaintextSpace.cpp CModulus.cpp FHEContext.cpp PAlgebra.cpp SingleCRT.cpp \
       DoubleCRT.cpp NumbTh.cpp bluestein.cpp IndexSet.cpp \
       Plaintext.cpp Util.cpp \
       FHE-SI.cpp Ciphertext.cpp \
       Serialization.cpp Matrix.cpp
       
TESTPROGS = Test_Regression_x Test_Statistics_x Test_General_x Test_AddMul_x

OBJPATHS = $(patsubst %.cpp,$(BUILD)/%.o, $(SRC))
TESTPATHS = $(addprefix $(TESTS)/, $(TESTPROGS))

SUFFIXES += .d
NODEPS = clean $(BUILD) $(TESTS) $(DISTR) obj
DEPFILES = $(patsubst %_x,$(TESTS)/%.d,$(TESTPROGS)) $(patsubst %.cpp,$(BUILD)/%.d,$(SRC))

all: $(OBJPATHS) $(TESTPATHS)

ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
    -include $(DEPFILES)
endif

obj: $(OBJPATHS)

test: $(OBJPATHS) $(TESTPATHS)

$(BUILD):
	mkdir -p $(BUILD)
$(TESTS):
	mkdir -p $(TESTS)

deps: $(DEPFILES)

$(BUILD)/%.d: %.cpp | $(BUILD)
	$(CC) $(CFLAGS) -MM -MT '$(BUILD)/$*.o' $< $(LDPATH) $(LDLIBS) -MF $@

$(TESTS)/%.d: %.cpp | $(TESTS)
	$(CC) $(CFLAGS) -MM -MT '$(TESTS)/$*_x' $< $(LDPATH) $(LDLIBS) -MF $@

$(BUILD)/%.o: %.cpp | $(BUILD)
	$(CC) $(CFLAGS) -o $@ -c $<

$(TESTS)/%_x: %.cpp $(OBJPATHS) $(TESTS)/%.d | $(TESTS)
	$(CC) $(CFLAGS) -o $@ $< $(LDPATH) $(OBJPATHS) $(LDLIBS)
  
clean:
	rm -rf $(BUILD) $(TESTS)/*_x $(TESTS)/*.d *~
