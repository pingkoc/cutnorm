from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()


setup(
    name='cutnorm',
    version='0.1.10',
    description=
    'Cutnorm approximation via Gaussian Rounding and Optimization with Orthogonality Constraints',
    long_description=readme(),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    url='https://github.com/pingkoc/cutnorm',
    author='Ping-Ko Chiu, Peter Diao, Olewasanmi Koyejo',
    author_email=
    'pchiu5@illinois.edu, peter.z.diao@gmail.com, sanmi@illinois.edu',
    license='MIT',
    packages=['cutnorm', 'cutnorm.tools'],
    install_requires=[
        'numpy>=1.16.0',
        'scikit-learn>=0.22.0',
        'scipy>=1.3.0',
        'numba>=0.45.0'
    ],
    python_requires='>=3.6',
    include_package_data=True,
    zip_safe=False)
