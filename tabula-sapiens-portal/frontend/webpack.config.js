const path = require("path");
const MiniCssExtractPlugin = require("mini-css-extract-plugin");

module.exports = {
  entry: {
    index: "./src/index.jsx",
  },
  output: {
    path: path.resolve("../backend/public/"),
    filename: "javascripts/[name].js",
  },
  resolve: {
    extensions: [".js", ".jsx"],
    alias: {
      "~": path.resolve(__dirname, "src"),
      styles: path.resolve(__dirname, "src/styles"),
    },
  },
  plugins: [
    new MiniCssExtractPlugin({
      filename: "styles/[name].bundle.min.css",
    }),
  ],
  module: {
    rules: [
      {
        test: /\.(js|jsx)$/,
        exclude: path.resolve(__dirname, "node_modules/"),
        use: {
          loader: "babel-loader",
        },
      },

      { test: /\.html$/, loader: "html-loader" },

      {
        test: /\.csv$/,
        loader: "csv-loader",
        options: {
          dynamicTyping: true,
          header: true,
          skipEmptyLines: true,
        },
      },
      {
        test: /\.(sa|sc|c)ss$/,
        exclude: path.resolve(__dirname, "node_modules/"),
        use: [
          MiniCssExtractPlugin.loader,
          {
            loader: "css-loader",
            options: {
              sourceMap: true,
              importLoaders: 2,
              modules: {
                mode: "local",
                localIdentName: "[local]-[hash:base64:5]",
              },
            },
          },
          {
            loader: "postcss-loader",
            options: {
              ident: "postcss",
              sourceMap: true,
              plugins: (loader) => [require("cssnano")({ preset: "default" })],
            },
          },
          {
            loader: "sass-loader",
            options: {
              sourceMap: true,
            },
          },
        ],
      },
      {
        test: /\.(png|svg|jpg|gif|pdf)$/,
        loader: "url-loader",
        options: {
          name: "[name].[ext]",
          publicPath: "assets/images/",
        },
      },
    ],
  },
  watch: true,
  optimization: {
    minimize: false,
  },
};
