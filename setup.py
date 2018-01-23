from setuptools import setup

setup(
    name='cutnorm',
    version='0.1.2',
    description=
    'Cutnorm approximation via Gaussian Rounding and Optimization with Orthogonality Constraints',
    url='https://github.com/pingkoc/cutnorm',
    author='Ping-Ko Chiu',
    author_email='pchiu5@illinois.edu',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    license='MIT',
    packages=['cutnorm'],
    install_requires=['numpy'],
    python_requires='>=3',
    zip_safe=False)
