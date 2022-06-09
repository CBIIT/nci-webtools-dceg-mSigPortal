const { Router } = require('express');

const router = Router();

router.get('/ping', (req, res) => res.send(true));

module.exports = router;
