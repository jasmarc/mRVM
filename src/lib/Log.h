// Copyright 2011 Jason Marcell

#ifndef SRC_LIB_LOG_H_
#define SRC_LIB_LOG_H_

#define VERBOSITY 2  // Set to 1, 2, 3 for norm, verbose, debug, 0 for silent

#define NORMAL     1
#define VERBOSE    2
#define DEBUG      3

#define LOG(L, format, ...) if (VERBOSITY >= L) fprintf (stderr, format, ## __VA_ARGS__);

#endif  // SRC_LIB_LOG_H_
