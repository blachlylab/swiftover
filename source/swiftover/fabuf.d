module swiftover.fabuf;

///
struct Fabuf(int windowSize)
{
    string contig;
    int start;
    int end;

    char[windowSize] buf;

    private void updateBuf();
    invariant
    {
        assert(this.end - this.start == windowSize);
    }
}