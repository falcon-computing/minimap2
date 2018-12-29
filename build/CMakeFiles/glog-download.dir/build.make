# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /curr/tianj/minimap2-ori/minimap2-mdbs

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /curr/tianj/minimap2-ori/minimap2-mdbs/build

# Utility rule file for glog-download.

# Include the progress variables for this target.
include CMakeFiles/glog-download.dir/progress.make

CMakeFiles/glog-download: CMakeFiles/glog-download-complete

CMakeFiles/glog-download-complete: deps/glog/src/glog-download-stamp/glog-download-install
CMakeFiles/glog-download-complete: deps/glog/src/glog-download-stamp/glog-download-mkdir
CMakeFiles/glog-download-complete: deps/glog/src/glog-download-stamp/glog-download-download
CMakeFiles/glog-download-complete: deps/glog/src/glog-download-stamp/glog-download-update
CMakeFiles/glog-download-complete: deps/glog/src/glog-download-stamp/glog-download-patch
CMakeFiles/glog-download-complete: deps/glog/src/glog-download-stamp/glog-download-configure
CMakeFiles/glog-download-complete: deps/glog/src/glog-download-stamp/glog-download-build
CMakeFiles/glog-download-complete: deps/glog/src/glog-download-stamp/glog-download-install
	$(CMAKE_COMMAND) -E cmake_progress_report /curr/tianj/minimap2-ori/minimap2-mdbs/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Completed 'glog-download'"
	/usr/bin/cmake -E make_directory /curr/tianj/minimap2-ori/minimap2-mdbs/build/CMakeFiles
	/usr/bin/cmake -E touch /curr/tianj/minimap2-ori/minimap2-mdbs/build/CMakeFiles/glog-download-complete
	/usr/bin/cmake -E touch /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-stamp/glog-download-done

deps/glog/src/glog-download-stamp/glog-download-install: deps/glog/src/glog-download-stamp/glog-download-build
	$(CMAKE_COMMAND) -E cmake_progress_report /curr/tianj/minimap2-ori/minimap2-mdbs/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "No install step for 'glog-download'"
	cd /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-build && /usr/bin/cmake -E touch /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-stamp/glog-download-install

deps/glog/src/glog-download-stamp/glog-download-mkdir:
	$(CMAKE_COMMAND) -E cmake_progress_report /curr/tianj/minimap2-ori/minimap2-mdbs/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Creating directories for 'glog-download'"
	/usr/bin/cmake -E make_directory /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/install
	/usr/bin/cmake -E make_directory /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-build
	/usr/bin/cmake -E make_directory /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog
	/usr/bin/cmake -E make_directory /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/tmp
	/usr/bin/cmake -E make_directory /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-stamp
	/usr/bin/cmake -E make_directory /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src
	/usr/bin/cmake -E touch /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-stamp/glog-download-mkdir

deps/glog/src/glog-download-stamp/glog-download-download: deps/glog/src/glog-download-stamp/glog-download-urlinfo.txt
deps/glog/src/glog-download-stamp/glog-download-download: deps/glog/src/glog-download-stamp/glog-download-mkdir
	$(CMAKE_COMMAND) -E cmake_progress_report /curr/tianj/minimap2-ori/minimap2-mdbs/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Performing download step (download, verify and extract) for 'glog-download'"
	cd /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog && /usr/bin/cmake -P /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-stamp/download-glog-download.cmake
	cd /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog && /usr/bin/cmake -P /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-stamp/verify-glog-download.cmake
	cd /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog && /usr/bin/cmake -P /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-stamp/extract-glog-download.cmake
	cd /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog && /usr/bin/cmake -E touch /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-stamp/glog-download-download

deps/glog/src/glog-download-stamp/glog-download-update: deps/glog/src/glog-download-stamp/glog-download-download
	$(CMAKE_COMMAND) -E cmake_progress_report /curr/tianj/minimap2-ori/minimap2-mdbs/build/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "No update step for 'glog-download'"
	/usr/bin/cmake -E touch /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-stamp/glog-download-update

deps/glog/src/glog-download-stamp/glog-download-patch: deps/glog/src/glog-download-stamp/glog-download-download
	$(CMAKE_COMMAND) -E cmake_progress_report /curr/tianj/minimap2-ori/minimap2-mdbs/build/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "No patch step for 'glog-download'"
	/usr/bin/cmake -E touch /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-stamp/glog-download-patch

deps/glog/src/glog-download-stamp/glog-download-configure: deps/glog/tmp/glog-download-cfgcmd.txt
deps/glog/src/glog-download-stamp/glog-download-configure: deps/glog/src/glog-download-stamp/glog-download-update
deps/glog/src/glog-download-stamp/glog-download-configure: deps/glog/src/glog-download-stamp/glog-download-patch
	$(CMAKE_COMMAND) -E cmake_progress_report /curr/tianj/minimap2-ori/minimap2-mdbs/build/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "No configure step for 'glog-download'"
	cd /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-build && /usr/bin/cmake -E touch /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-stamp/glog-download-configure

deps/glog/src/glog-download-stamp/glog-download-build: deps/glog/src/glog-download-stamp/glog-download-configure
	$(CMAKE_COMMAND) -E cmake_progress_report /curr/tianj/minimap2-ori/minimap2-mdbs/build/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "No build step for 'glog-download'"
	cd /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-build && /usr/bin/cmake -E touch /curr/tianj/minimap2-ori/minimap2-mdbs/build/deps/glog/src/glog-download-stamp/glog-download-build

glog-download: CMakeFiles/glog-download
glog-download: CMakeFiles/glog-download-complete
glog-download: deps/glog/src/glog-download-stamp/glog-download-install
glog-download: deps/glog/src/glog-download-stamp/glog-download-mkdir
glog-download: deps/glog/src/glog-download-stamp/glog-download-download
glog-download: deps/glog/src/glog-download-stamp/glog-download-update
glog-download: deps/glog/src/glog-download-stamp/glog-download-patch
glog-download: deps/glog/src/glog-download-stamp/glog-download-configure
glog-download: deps/glog/src/glog-download-stamp/glog-download-build
glog-download: CMakeFiles/glog-download.dir/build.make
.PHONY : glog-download

# Rule to build all files generated by this target.
CMakeFiles/glog-download.dir/build: glog-download
.PHONY : CMakeFiles/glog-download.dir/build

CMakeFiles/glog-download.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/glog-download.dir/cmake_clean.cmake
.PHONY : CMakeFiles/glog-download.dir/clean

CMakeFiles/glog-download.dir/depend:
	cd /curr/tianj/minimap2-ori/minimap2-mdbs/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /curr/tianj/minimap2-ori/minimap2-mdbs /curr/tianj/minimap2-ori/minimap2-mdbs /curr/tianj/minimap2-ori/minimap2-mdbs/build /curr/tianj/minimap2-ori/minimap2-mdbs/build /curr/tianj/minimap2-ori/minimap2-mdbs/build/CMakeFiles/glog-download.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/glog-download.dir/depend

