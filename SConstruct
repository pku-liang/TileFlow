import os 


AddOption('--static', dest='link_static', default=False, action='store_true', help='Use static linking (default is dynamic)')
AddOption('--accelergy', dest='use_accelergy', default=False, action='store_true', help='Build Timeloop with Accelergy (default is to use pat/src)')
AddOption('--d', dest='debug', default=False, action='store_true', help='Debug build (default is off)')
AddOption('--with-isl', dest='with_isl', default=False, action='store_true', help='Build with ISL support (default is false)')
AddOption('--clang', dest='clang', default=(str(Platform())=='darwin'), action='store_true', help='Build using clang (default is true for MacOS, otherwise false)')
AddOption('--parser', dest='parser', default=False, action = 'store_true', help='Build frontend parser')

if GetOption("parser"):
  print ('parser build') 
  VariantDir('build', 'tests', duplicate=0)
else: 
  print ('application build') 
  VariantDir('build', 'src', duplicate = 0)

env = Environment(ENV = os.environ)


if GetOption('clang'):
  print('Building with clang instead of gcc.')

if not GetOption('clang'):
  env.Replace(AR = "gcc-ar")
  env.Replace(RANLIB = "gcc-ranlib")

env.Append(TIMELOOP_BASE_DIR = Dir('./3rdparty/timeloop/').abspath)
env.Append(BUILD_BASE_DIR = Dir('.').abspath)
env.SConscript('build/SConscript', exports='env')
