#include <iostream>
#include <pybind11/embed.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace py::literals;

int main() {
    // start the interpreter and keep it alive
    py::scoped_interpreter guard{};

    // use the Python API
    py::print("Hello, World!");

    // executing code from string
    py::exec(R"(
        kwargs = dict(name="World", number=42)
        message = "Hello, {name}! The answer is {number}".format(**kwargs)
        print(message)
        my_variable = 20
    )");

    // evaluating expressions in python
    // (you can also use py::eval_file("script.py");)
    int result = py::eval("my_variable + 10").cast<int>();
    py::print(result);

    // use python with the API
    auto locals = py::dict("name"_a="World", "number"_a=42);
    py::exec(R"(
        message = "Hello, {name}! The answer is {number}".format(**locals())
    )", py::globals(), locals);
    auto message = locals["message"].cast<std::string>();
    py::print(message);

    // importing modules
    py::module sys = py::module::import("sys");
    py::print(sys.attr("path"));

    // importing local file as module
    // (the add function is seen as an attribute of the module)
    py::module calc = py::module::import("calc");
    py::object add_function = calc.attr("add");
    py::object add_result = add_function(1, 2);
    int n = add_result.cast<int>();
    assert(n == 3);
    py::print(add_result.cast<int>());

    // plotting with python
    // https://matplotlib.org/gallery/lines_bars_and_markers/simple_plot.html#sphx-glr-gallery-lines-bars-and-markers-simple-plot-py
    // import modules
    py::module matplotlib = py::module::import("matplotlib");
    py::module plt = py::module::import("matplotlib.pyplot");
    py::module np = py::module::import("numpy");
    std::cout << py::str(matplotlib) << std::endl;
    // failing to find a module
    try {
        py::module fail_module = py::module::import("fake_module_name_that_doesnt_exist");
    } catch (...) {
        std::cout << "You can't import any name you want" << std::endl;
    }
    // create data we want to plot
    // t = [0.00, 0.01, 0.02, ..., 2.00]
    std::vector<double> t(static_cast<size_t>(2.0/0.01));
    std::generate(t.begin(), t.end(), [range = 0.0] () mutable { return range += 0.01; });
    // s = 1 + sin(pi * 2 * t)
    std::vector<double> s(t.size());
    // get pi from python
    // (just as an example - of course we could use C++ for that)
    auto pi = np.attr("pi").cast<double>();
    // get the sin function from python
    // (just as an example - of course we could use C++ for that
    auto sin = [&np](double x) {
        static py::object py_sin = np.attr("sin");
        return py_sin(x).cast<double>();
    };
    // s = 1 + sin(pi * 2 * t)
    std::transform(t.begin(),t.end(),s.begin(),[&sin,&pi](double x){return 1 + sin(pi * 2 * x);});
    // use the plotting library (subplot returns a tuple of return values)
    py::tuple pytuple = plt.attr("subplots")();
    py::object fig = pytuple[0];
    py::object ax = pytuple[1];
    // interactive mode
    plt.attr("ion")();
    // plot t and s (pybind automatically converts t and s to python types)
    ax.attr("plot")(t, s);
    // set labels
    ax.attr("set")("xlabel"_a="time (s)", "ylabel"_a="voltage (mV)", "title"_a="About as simple as it gets, folks");
    // put grid behind plot
    ax.attr("grid")();
    // save a file with the figure
    fig.attr("savefig")("test.png");
    // draw on screen
    plt.attr("draw")();
    plt.attr("pause")(0.0001);
    for (int i = 0; i < s.size(); ++i) {
        // change values and plot again a few times
        std::rotate(s.begin(),s.begin()+1,s.end());
        plt.attr("cla")();
        ax.attr("plot")(t, s);
        ax.attr("set")("xlabel"_a="time (s)", "ylabel"_a="voltage (mV)", "title"_a="About as simple as it gets, folks");
        plt.attr("draw")();
        plt.attr("pause")(0.0001);
    }
    plt.attr("show")();
    return 0;
}