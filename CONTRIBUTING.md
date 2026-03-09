# Contributing to EarlGrey ParTEA

We welcome contributions to EarlGrey ParTEA!

## How to Contribute

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add tests if applicable
5. Commit your changes (`git commit -m 'Add amazing feature'`)
6. Push to the branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

## Reporting Issues

Please use the [GitHub issue tracker](https://github.com/TobyBaril/EarlGreyParTEA/issues) to:
- Report bugs
- Request features
- Ask questions

## Development Setup

```bash
git clone https://github.com/TobyBaril/EarlGreyParTEA.git
cd EarlGreyParTEA
chmod +x earlGreyParTEA*
export PATH="$PWD:$PATH"
```

## Code Style

- Follow PEP 8 for Python code
- Use descriptive variable names
- Add comments for complex logic
- Update documentation for new features

## Testing

Before submitting a PR, please test your changes:

```bash
# Test wrapper scripts
./earlGreyParTEA --help
./earlGreyParTEA --generate-config test.yaml

# Test dry run
earlGreyParTEA -c test.yaml -t 4 --dry-run
```

## Questions?

Contact: tobias.baril[at]unine.ch

Thank you for contributing!
