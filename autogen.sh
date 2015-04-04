#!/bin/sh -e
test -n "$srcdir" || srcdir=`dirname "$0"`
test -n "$srcdir" || srcdir=.

autoreconf --force --install --verbose "$srcdir"
aclocal --install -I m4
automake --add-missing
test -n "$NOCONFIGURE" || "$srcdir/configure" "$@"
