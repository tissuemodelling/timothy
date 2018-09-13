#* **************************************************************************
# * This file is part of Timothy
# *
# * Copyright (c) 2014/15 Maciej Cytowski
# * Copyright (c) 2014/15 ICM, University of Warsaw, Poland
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# *
# * *************************************************************************/

.PHONY: all clean 

all:
	@ if command -v gmake >/dev/null 2>&1; then \
	echo "using gmake to build timothy"; \
	cd src/; gmake; \
        else \
        echo "using make to build timothy"; \
        cd src/; make; \
        fi

install:
	@ if command -v gmake >/dev/null 2>&1; then \
	echo "installing timothy"; \
	cd src/; gmake install; \
	else \
	echo "installing timothy"; \
	cd src/; make install; \
	fi

clean:
	@ if command -v gmake >/dev/null 2>&1; then \
        echo "using gmake to clean timothy"; \
        cd src/; gmake clean; \
        else \
        echo "using make to clean timothy"; \
        cd src/; make clean; \
        fi
