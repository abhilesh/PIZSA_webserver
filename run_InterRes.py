#!/usr/bin/env python

"""
Copyright (C) 2016 Abhilesh Dhawanjewar <abhilesh7@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Abhilesh Dhawanjewar and Neelesh Soni"
__copyright__ = "Copyright C 2016"
__license__ = "LGPL"
__version__ = "1.0"
__maintainer__ = "Abhilesh Dhawanjewar"
__email__ = "abhilesh7@gmail.com"

import scripts.src_InterRes.main_InterRes
import time

start_time = time.time()
scripts.src_InterRes.main_InterRes.main_InterRes()
end_time = time.time()

print "Time elapsed: ", round(end_time - start_time, 3), 's'
