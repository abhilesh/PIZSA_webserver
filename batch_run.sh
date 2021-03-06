#!/bin/bash

# Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com>

# This file is part of PIZSA.

# PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.


for i in `ls /home/abhilesh/train_set/ | grep .pdb$`;
do
    echo $i" Running";
	./run_PIZSA.py /home/abhilesh/train_set/$i;
done;
