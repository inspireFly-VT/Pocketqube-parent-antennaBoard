SRC=main
TARGET=patch-doc-v1.0

ifndef BUILD_DIR
	BUILD_DIR=build
endif

PDF_READER=evince

all:
	mkdir -p $(BUILD_DIR)
	mkdir -p $(BUILD_DIR)/references
	cp -u references/*.bib $(BUILD_DIR)/references/ 2>/dev/null || true
	latexmk -pdf -jobname=$(BUILD_DIR)/$(TARGET) $(SRC).tex
	-makeindex $(BUILD_DIR)/$(TARGET).nlo -s nomencl.ist -o $(BUILD_DIR)/$(TARGET).nls 2>/dev/null || true
	latexmk -pdf -jobname=$(BUILD_DIR)/$(TARGET) $(SRC).tex
	rm -f $(TARGET).pdf
	cp $(BUILD_DIR)/$(TARGET).pdf .

index:
	makeindex $(BUILD_DIR)/$(TARGET).nlo -s nomencl.ist -o $(BUILD_DIR)/$(TARGET).nls

view:
	$(PDF_READER) $(BUILD_DIR)/$(TARGET).pdf &

clean:
	rm $(BUILD_DIR)/* header/*.aux chapters/*.aux references/*.aux
