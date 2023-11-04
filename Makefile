# Targets
.PHONY: all build run clean distclean 

# Default target (all): Run all the main tasks
all: build

# Target: build
# Build the functions and directories needed
build:
	mkdir -p images results data
	
# Target: run
# Run tests
run: build
  
# Target: clean
# Clean intermediate files
clean:
	@$(RM) *.aux *.log *.pdf *.txt
	@$(RM) .Rhistory
	
# Target: distclean
# Clean all generated files and directories
distclean: clean
	@$(RM) data/ -r
	@$(RM) images/ -r
	@$(RM) results/ -r

