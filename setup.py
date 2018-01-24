from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()


setup(
    name='cutnorm',
    version='0.1.5',
    description=
    'Cutnorm approximation via Gaussian Rounding and Optimization with Orthogonality Constraints',
    long_description=readme(),
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
    url='https://github.com/pingkoc/cutnorm',
    author='Ping-Ko Chiu',
    author_email='pchiu5@illinois.edu',
    license='MIT',
    packages=['cutnorm', 'cutnorm.tools'],
    install_requires=['numpy'],
    python_requires='>=3',
    include_package_data=True,
    zip_safe=False)
