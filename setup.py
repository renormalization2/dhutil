from setuptools import setup, find_packages

setup(
    name="dhutil",
    version="0.1",
    packages=find_packages(),  # Automatically find and include all packages
    install_requires=[  # List of dependencies
        "numpy",
        "astropy",
        # "matplotlib",
        # 'requests', etc.
    ],
    entry_points={  # Optional: entry points for command-line scripts
        "console_scripts": [
            # 'my_command = my_project.module1:main_function',  # Replace with actual command and function
        ],
    },
    author="dhhyun",
    author_email="hdhd333@gmail.com",
    description="dhutil",
    url="https://github.com/username/my_project",  # Replace with your GitHub URL
    classifiers=[
        "Programming Language :: Python :: 3",
        # "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
