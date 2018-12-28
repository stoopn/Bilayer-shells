# $Id: Makefile 2501 2007-11-20 02:33:29Z benkirk $


# The location of the mesh library
meshdir := . 

# include the library options determined by configure.  This will
# set the variables INCLUDE and LIBS that we will need to build and
# link with the library.
include /home/libs/LIBMESHSUBDIV/Make.common
 
#
###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# source files
srcfiles       := shellelement.C shellsystem.C ExplicitShellsystem.C Tri_Tri_Intersection.C VoxelList.C shellmeasurementsystem.C

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
expl_obj	:= $(patsubst %.C, %.$(obj-suffix), explicit_shell.C)

###############################################################################



.PHONY: clean clobber distclean

###############################################################################
# Target:
#
explicit 	   := ./explicitshell-$(METHOD)
# steady		   := ./steadyshell-$(METHOD)

expl:: $(explicit)

$(expl_obj): explicit_shell.C
	@echo "Compiling explicit_shell.C..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(libmesh_INCLUDE) -c explicit_shell.C  -o $(expl_obj) 

$(explicit): $(objects) $(expl_obj)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects) $(expl_obj) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS) -lvtkFiltering  

# Useful rules.
clean:
	@rm -f $(objects) *~

clobber:
	@$(MAKE) clean
	@rm -f $(target) out.gmv

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o

# Warning, the 3D problem may be extremely slow if you are running in debug mode.
run: $(target)
	@echo "***************************************************************"
	@echo "* Running Example " $(LIBMESHRUN) $(target) $(LIBMESHOPTIONS)
	@echo "***************************************************************"
	@echo " "
	@$(LIBMESHRUN) $(target) -d 1 -n 20 $(LIBMESHOPTIONS)
	@$(LIBMESHRUN) $(target) -d 2 -n 15 $(LIBMESHOPTIONS)
	@$(LIBMESHRUN) $(target) -d 3 -n 6 $(LIBMESHOPTIONS)
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Example " $(LIBMESHRUN) $(target) $(LIBMESHOPTIONS)
	@echo "***************************************************************"


# include the dependency list
include .depend


#
# Dependencies
#
.depend:
	@$(perl) /home/libs/LIBMESHSUBDIV/contrib/bin/make_dependencies.pl -I. $(foreach i, $(wildcard $(meshdir)/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend
	@echo "Updated .depend"

###############################################################################
