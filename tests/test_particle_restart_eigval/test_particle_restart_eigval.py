#!/usr/bin/env python

import sys
sys.path.insert(0, '..')
from testing_harness import ParticleRestartTestHarness


if __name__ == '__main__':
    harness = ParticleRestartTestHarness('particle_9_555.*')
    harness.main()
