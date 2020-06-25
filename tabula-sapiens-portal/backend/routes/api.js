var express = require("express");
var router = express.Router();

var app = express();
var multer = require("multer");
var cors = require("cors");

/* Upload file */
router.get("/*", function (req, res, next) {
  app.use(cors());
  var storage = multer.diskStorage({
    destination: function (req, file, cb) {
      cb(null, "public/annotations");
    },
    filename: function (req, file, cb) {
      cb(null, Date.now() + "-" + file.originalname);
    },
  });
  app.post("/upload", function (req, res) {
    upload(req, res, function (err) {
      if (err instanceof multer.MulterError) {
        return res.status(500).json(err);
      } else if (err) {
        return res.status(500).json(err);
      }
      return res.status(200).send(req.file);
    });
  });
});

module.exports = router;
