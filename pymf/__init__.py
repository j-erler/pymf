# -*- coding: utf-8 -*-

# 
#  This file is part of pyMF.
# 
#  pyMF is free software; you can redistribute it and/or modify
#  it under the terms of the MIT License.
# 
#  pyMF is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  provided copy of the MIT License for more details.
# 
# 
"""pyMF is a python implementation of matched filtering techniques for 
   applications in astronomy.
"""

__version__ = "1.0"

__bibtex__ = """

"""

from .pymf import (corr2cov, rad_profile, power_spec, cross_spec, make_filter_map, filter_map_mf, filter_map_cmf, filter_map_mmf, filter_map_cmmf)
