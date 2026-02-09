# Contributing Guidelines

Thanks for taking the time to contribute to BloodFlowTrixi.jl! ❤️

We appreciate your interest in improving this project. Before you start contributing, please take a moment to review the following guidelines.


## Reporting Issues

First off, we assume that you have read the available [Documentation](https://juliahealth.org/BloodFlowTrixi.jl/dev/).

Before you report an issue, it is best to search for existing [Issues](https://github.com/JuliaHealth/BloodFlowTrixi.jl/issues) that might help you. In case you have found a suitable issue and still need clarification, you can write your question in this issue.

If you then still feel the need to report an issue, we recommend the following:

- Open an [Issue on GitHub](https://github.com/JuliaHealth/BloodFlowTrixi.jl/issues/new).
- Provide as much context as you can about what you're running into.
- Provide project and platform versions, depending on what seems relevant.

We will then take care of the issue as soon as possible.


## How to Contribute

- **Implement Code Changes**: Introduce your modifications, ensuring adherence to the [Julia Blue Style Guidelines](https://github.com/invenia/BlueStyle). For new features, include informative comments, docstrings, and consider enriching the documentation with relevant examples.

- **Format Code**: This project uses JuliaFormatter with Blue style. Run formatting from the repository root using these steps:
	- `julia --project=dev`
	- `using JuliaFormatter`
	- `format(".")`

- **Validate Changes with Tests**: Execute existing tests to verify the compatibility of your alterations with the current functionality. If applicable, incorporate additional tests to validate your new contributions.

- **Undergo Code Review**: Subject your code to review by maintainers, who will provide feedback and may request further adjustments before merging.


## License

By contributing to BloodFlowTrixi.jl, you agree that your contributions will be licensed under the [MIT License](LICENSE).

Thank you for contributing to BloodFlowTrixi.jl! 🌟
