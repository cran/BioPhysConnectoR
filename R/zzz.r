#############################################
#   This code is subject to the license as stated in DESCRIPTION.
#   Using this software implies acceptance of the license terms:
#    - GPL 2
#
#   (C) by F. Hoffgaard, P. Weil, and K. Hamacher in 2009.
#
#  hoffgaard(AT)bio.tu-darmstadt.de
#
#
#  http://www.kay-hamacher.de
#############################################



.First.lib <- function(lib, pkg) {
  library.dynam("BioPhysConnectoR", pkg, lib)
}
